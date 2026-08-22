import numpy as np
import matplotlib.pyplot as plt
import os, sys, re

# ============================================================
# USER SETTINGS
# ============================================================
N_BINS    = 8
REMOVE_N  = 2
COMPONENT = 'JPz'   # JPx | JPy | JPz  (transport direction!)

IFACE_LEFT_0IDX  = 3   # group 4 — last Bi2Te3 bin (0-indexed)
IFACE_RIGHT_0IDX = 4   # group 5 — first Sb2Te3 bin (0-indexed)
IFACE_LABEL      = 'Bi₂Te₃ | Sb₂Te₃'

# None = ALL rows | e.g. 2500 = last 2500 rows only
JP_STEPS  = None
# ============================================================
dt_fs           = 1.0
sample_interval = 10
output_interval = 100
time_per_row_ns = dt_fs * sample_interval * output_interval * 1e-6
eV_to_J         = 1.602176634e-19
amu_to_kg        = 1.66053906660e-27
JP_conv          = (eV_to_J ** 1.5) * (amu_to_kg ** -0.5)

FIT_0IDX = list(range(REMOVE_N, N_BINS - REMOVE_N))   # [2,3,4,5]
FIT_GIDS = [i + 1 for i in FIT_0IDX]                  # [3,4,5,6]

LEFT_FIT_0IDX  = list(range(REMOVE_N, IFACE_LEFT_0IDX + 1))   # [2,3]
RIGHT_FIT_0IDX = list(range(IFACE_RIGHT_0IDX, N_BINS - REMOVE_N))  # [4,5]

# ============================================================
# HELPERS
# ============================================================

def parse_relaxed_xyz(folder):
    fpath = os.path.join(folder, 'relaxed.xyz')
    if not os.path.exists(fpath):
        print(f"  WARNING: relaxed.xyz not found, falling back to nemd_setup.txt")
        return None, None, {}

    with open(fpath, encoding='utf-8', errors='replace') as f:
        lines = f.readlines()

    n_atoms = int(lines[0].strip())
    header  = lines[1]

    lm = re.search(r'Lattice="([^"]+)"', header)
    if not lm:
        print("  WARNING: no Lattice= in relaxed.xyz header")
        return None, None, {}
    lvals = list(map(float, lm.group(1).split()))
    cell  = np.array(lvals).reshape(3, 3)

    Lz_ang = cell[2, 2]
    A_ang2 = float(np.linalg.norm(np.cross(cell[0], cell[1])))

    has_group = 'group' in header
    syms, pos, groups = [], [], []
    for line in lines[2:2 + n_atoms]:
        parts = line.split()
        if len(parts) < 4:
            continue
        syms.append(parts[0])
        pos.append([float(parts[1]), float(parts[2]), float(parts[3])])
        groups.append(int(parts[4]) if (has_group and len(parts) > 4) else -1)

    pos    = np.array(pos)
    groups = np.array(groups)

    centers_nm = {}
    for g in range(1, N_BINS + 1):
        mask = groups == g
        if mask.sum() == 0:
            continue
        z_min = pos[mask, 2].min()
        z_max = pos[mask, 2].max()
        centers_nm[g] = (z_min + z_max) / 2.0 * 0.1  # Ang -> nm

    print(f"  relaxed.xyz: Lz={Lz_ang:.4f} Ang  A={A_ang2:.4f} Ang^2")
    print(f"  Bin centers from relaxed.xyz (nm): {centers_nm}")
    return Lz_ang, A_ang2, centers_nm


def parse_setup_fallback(folder):
    L, A = None, None
    try:
        with open(os.path.join(folder, 'nemd_setup.txt'), encoding='utf-8', errors='replace') as f:
            for line in f:
                s = line.strip()
                if L is None and re.search(r'(Lz\s*\(total\)|L_transport)', s) and ':' in s:
                    m = re.search(r':\s*([\d.]+)', s)
                    if m: L = float(m.group(1))
                if A is None and 'A_cross' in s and ':' in s:
                    m = re.search(r':\s*([\d.]+)', s)
                    if m: A = float(m.group(1))
    except Exception as e:
        print(f"  WARNING parse_setup_fallback: {e}")
    return L, A


def read_bin_centers_fallback(folder):
    centers = {}
    try:
        with open(os.path.join(folder, 'nemd_setup.txt'), encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
        in_table = False
        gid_seen = {}
        for line in lines:
            s = line.strip()
            if 'GROUP BOUNDARIES' in s:
                in_table = True; continue
            if in_table:
                if s.startswith('===') or 'L_eff' in s:
                    break
                if s.startswith('GID') or re.match(r'^-+$', s):
                    continue
                s_clean = re.split(r'←|<-|#', s)[0]
                nums = re.findall(r'[\d]+\.[\d]+|[\d]+', s_clean)
                if len(nums) >= 5:
                    gid   = int(nums[0])
                    z_min = float(nums[-5])
                    z_max = float(nums[-4])
                    if 1 <= gid <= N_BINS and gid not in gid_seen:
                        centers[gid] = (z_min + z_max) / 2.0 * 0.1
                        gid_seen[gid] = True
    except Exception as e:
        print(f"  WARNING read_bin_centers_fallback: {e}")
    return centers


def parse_reservoir_bounds(folder):
    """Read the ACTUAL physical extent (in nm) of the thermostatted source
    and sink regions from the GROUP BOUNDARIES table."""
    src_bounds, snk_bounds = None, None
    try:
        with open(os.path.join(folder, 'nemd_setup.txt'), encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
        for line in lines:
            tok = line.split()
            if len(tok) < 5:
                continue
            role = tok[1].lower() if len(tok) > 1 else ''
            if role == 'src':
                src_bounds = (float(tok[3]) * 0.1, float(tok[4]) * 0.1)  # Ang -> nm
            elif role == 'snk':
                snk_bounds = (float(tok[3]) * 0.1, float(tok[4]) * 0.1)
    except Exception as e:
        print(f"  WARNING parse_reservoir_bounds: {e}")
    return src_bounds, snk_bounds


def fit_gradient(x_m, T_vals):
    coeffs = np.polyfit(x_m, T_vals, 1)
    T_pred = np.polyval(coeffs, x_m)
    ss_res = np.sum((T_vals - T_pred) ** 2)
    ss_tot = np.sum((T_vals - T_vals.mean()) ** 2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return coeffs[0], coeffs[1], R2


def get_jp_col_start(ncols):
    comp      = COMPONENT.lower()
    block_idx = {'jpx': 0, 'jpy': 1, 'jpz': 2, 'jkx': 3, 'jky': 4, 'jkz': 5}[comp]
    T_end     = 1 + N_BINS
    jp_start  = T_end + block_idx * (N_BINS + 1) + 1
    return jp_start


# ============================================================
# COMPUTE ONE FOLDER
# ============================================================
def compute_folder(folder):
    fpath = os.path.join(folder, 'compute.out')
    if not os.path.exists(fpath):
        print(f"  SKIP {folder} — compute.out not found")
        return None

    L_ang, A_ang2, centers_nm = parse_relaxed_xyz(folder)

    if L_ang is None or A_ang2 is None:
        print("  Falling back to nemd_setup.txt for L and A...")
        L_ang, A_ang2 = parse_setup_fallback(folder)

    if not centers_nm:
        print("  Falling back to nemd_setup.txt for bin centers...")
        centers_nm = read_bin_centers_fallback(folder)

    src_bounds_nm, snk_bounds_nm = parse_reservoir_bounds(folder)
    if src_bounds_nm and snk_bounds_nm:
        print(f"  Reservoirs (nm): src {src_bounds_nm[0]:.3f}-{src_bounds_nm[1]:.3f}  "
              f"snk {snk_bounds_nm[0]:.3f}-{snk_bounds_nm[1]:.3f}")

    if L_ang is None or A_ang2 is None:
        print(f"  SKIP {folder} — geometry parse failed")
        return None

    print(f"  Using geometry: L={L_ang:.4f} Ang  A={A_ang2:.4f} Ang^2")

    data         = np.loadtxt(fpath)
    ndata, ncols = data.shape
    t_ns         = np.arange(ndata) * time_per_row_ns

    jp_start_row = 0 if JP_STEPS is None else max(0, ndata - JP_STEPS)
    steps_jp     = ndata - jp_start_row

    A_m2  = A_ang2 * 1e-20
    L_nm  = L_ang  * 0.1
    A_nm2 = A_ang2 * 1e-2

    cvals = [centers_nm.get(g, np.nan) for g in range(2, N_BINS)]
    cvals = [v for v in cvals if not np.isnan(v)]
    if len(cvals) >= 2:
        L_bin_nm  = float(np.mean(np.diff(cvals)))
        L_bin_ang = L_bin_nm * 10.0
        L_bin_m   = L_bin_ang * 1e-10
        print(f"  L_bin (mean spacing G2-G7): {L_bin_nm:.4f} nm")
    else:
        L_bin_ang = L_ang / N_BINS
        L_bin_m   = L_bin_ang * 1e-10
        L_bin_nm  = L_bin_ang * 0.1

    Vol_nm3 = A_nm2 * L_bin_nm

    if len(centers_nm) < N_BINS:
        print(f"  WARNING: only {len(centers_nm)}/{N_BINS} bin centers parsed")

    x_all_nm = np.array([centers_nm.get(g, np.nan) for g in range(1, N_BINS + 1)])
    x_all_m  = x_all_nm * 1e-9

    # Temperature
    T_all  = data[jp_start_row:, 1:N_BINS + 1]
    T_mean = np.mean(T_all, axis=0)

    print(f"  T_mean per bin: " +
          "  ".join([f"G{i+1}:{T_mean[i]:.2f}K" for i in range(N_BINS)]))

    x_fit_m = x_all_m[FIT_0IDX]
    T_fit   = T_mean[FIT_0IDX]
    grad_T, intercept, R2_T = fit_gradient(x_fit_m, T_fit)
    deltaT = T_mean[0] - T_mean[N_BINS - 2]

    jp_col  = get_jp_col_start(ncols)
    JP_all  = data[jp_start_row:, jp_col:jp_col + N_BINS]
    JP_mean = np.mean(JP_all, axis=0)
    avg_jp  = np.mean(JP_mean[FIT_0IDX])
    Q_Wm2   = (avg_jp * JP_conv) / (A_m2 * L_bin_m)
    kappa   = -Q_Wm2 / grad_T

    x_left_m = x_all_m[LEFT_FIT_0IDX]
    T_left   = T_mean[LEFT_FIT_0IDX]
    grad_L, intc_L, R2_L = fit_gradient(x_left_m, T_left)

    x_right_m = x_all_m[RIGHT_FIT_0IDX]
    T_right   = T_mean[RIGHT_FIT_0IDX]
    grad_R, intc_R, R2_R = fit_gradient(x_right_m, T_right)

    x_iface_m     = 0.5 * (x_all_m[IFACE_LEFT_0IDX] + x_all_m[IFACE_RIGHT_0IDX])
    T_iface_L     = grad_L * x_iface_m + intc_L
    T_iface_R     = grad_R * x_iface_m + intc_R
    delta_T_iface = T_iface_L - T_iface_R

    avg_jp_rk = np.mean(JP_mean[FIT_0IDX])
    Q_rk      = (avg_jp_rk * JP_conv) / (A_m2 * L_bin_m)
    R_K       = delta_T_iface / abs(Q_rk)
    R_K_nm2KW = R_K * 1e9
    G_K       = 1.0 / R_K if R_K != 0 else np.inf

    label = os.path.basename(os.path.abspath(folder)) if folder == '.' else folder
    print(f"\n  {label}  |  L={L_nm:.2f}nm  A={A_nm2:.3f}nm²")
    print(f"  Q = {Q_Wm2:.4e} W/m²  ({COMPONENT})")
    print(f"  Global:  dT/dx={grad_T:.4e} K/m  kappa={kappa:.4f} W/mK  R2={R2_T:.4f}")
    print(f"  Left  (Bi2Te3):  dT/dx={grad_L:.4e}  T@iface={T_iface_L:.4f} K  R2={R2_L:.4f}")
    print(f"  Right (Sb2Te3):  dT/dx={grad_R:.4e}  T@iface={T_iface_R:.4f} K  R2={R2_R:.4f}")
    print(f"  dT_interface = {delta_T_iface:.4f} K")
    print(f"  R_K = {R_K:.4e} m2K/W  =  {R_K_nm2KW:.4f} x1e-9 m2K/W")
    print(f"  G_K = {G_K:.4e} W/m2K  =  {G_K/1e6:.4f} MW/m2K")

    return dict(
        folder=label, L_nm=L_nm, A_nm2=A_nm2, A_m2=A_m2,
        L_bin_nm=L_bin_nm, Vol_nm3=Vol_nm3,
        ndata=ndata, steps_jp=steps_jp,
        T_all=T_all,
        T_mean=T_mean, T_fit=T_fit,
        x_all_nm=x_all_nm, x_all_m=x_all_m,
        x_fit_m=x_fit_m, centers_nm=centers_nm,
        src_bounds_nm=src_bounds_nm, snk_bounds_nm=snk_bounds_nm,
        grad_T=grad_T, intercept=intercept, R2_T=R2_T, deltaT=deltaT,
        grad_L=grad_L, intc_L=intc_L, R2_L=R2_L,
        grad_R=grad_R, intc_R=intc_R, R2_R=R2_R,
        x_iface_m=x_iface_m,
        T_iface_L=T_iface_L, T_iface_R=T_iface_R,
        delta_T_iface=delta_T_iface,
        JP_mean=JP_mean, JP_all=JP_all, avg_jp=avg_jp,
        Q_Wm2=Q_Wm2, Q_rk=Q_rk, kappa=kappa,
        R_K=R_K, R_K_nm2KW=R_K_nm2KW, G_K=G_K,
        t_ns=t_ns,
    )


# ============================================================
# JOURNAL-QUALITY GLOBAL SETTINGS
# ============================================================
plt.rcParams.update({
    "font.family":       "serif",
    "font.serif":        ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset":  "stix",
    "font.size":         16,
    "axes.labelsize":    18,
    "axes.titlesize":    18,
    "xtick.labelsize":   18,
    "ytick.labelsize":   18,
    "legend.fontsize":   14,
    "axes.linewidth":    2.0,
    "lines.linewidth":   2.5,
    "lines.markersize":  8,
    "xtick.major.width": 2.0,
    "ytick.major.width": 2.0,
    "xtick.minor.width": 1.5,
    "ytick.minor.width": 1.5,
})

def apply_style(ax):
    ax.tick_params(axis='both', which='major',
                   direction='inout', length=8, width=2.0,
                   top=False, right=False)
    ax.tick_params(axis='both', which='minor',
                   direction='in', length=6, width=1.5,
                   top=False, right=False)
    ax.minorticks_on()


# ============================================================
# JP vs BIN PLOT  (same style as post_processing_nemd_Jp.py)
# ============================================================
def plot_jp_vs_bin(common, results):
    """
    One subplot per folder — all bins as solid blue line+dots,
    fit bins as open red circles, dashed red avg line.
    Title: folder name + Q (W/m²).  Saved as JP_vs_bin_<COMPONENT>.pdf/png
    """
    n = len(common)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4))
    if n == 1:
        axes = [axes]

    for ax, f in zip(axes, common):
        v    = results[f]
        bins = np.arange(1, N_BINS + 1)

        # plot -JP so positive = heat flowing in transport direction
        ax.plot(bins, -v['JP_mean'], 'o-', color='C0', lw=2, ms=8, label='All bins')
        ax.plot(
            [bins[i] for i in FIT_0IDX],
            -v['JP_mean'][FIT_0IDX],
            'ro', ms=10, mfc='white', mew=2, label='Fit bins'
        )
        ax.axhline(
            -v['avg_jp'], color='C3', ls='--', lw=1.5,
            label=f'avg={-v["avg_jp"]:.4f}'
        )
        ax.axhline(0, color='k', lw=0.6, ls=':')

        ax.set_xlabel('Bin index', fontsize=13)
        ax.set_ylabel(
            f'{COMPONENT} (eV$^{{3/2}}$ amu$^{{-1/2}}$)', fontsize=12
        )
        ax.set_title(f'{f}  Q={-v["Q_Wm2"]:.3e} W/m²', fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fname = f"JP_vs_bin_{COMPONENT}.{ext}"
        plt.savefig(fname, dpi=300, bbox_inches='tight')
        print(f"  Saved: {fname}")
    plt.show()
    plt.close()


# ============================================================
# JP DEVIATION FROM MEAN PLOT
# ============================================================
def plot_jp_deviation(common, results):
    """
    Two-panel layout per folder:
      Top: raw JP vs bin (same as plot_jp_vs_bin)
      Bottom: deviation of each bin's JP from the fit-bin mean,
              coloured by sign (red = above mean, blue = below),
              with ±1σ band across fit bins and per-bin std bars.

    Physical meaning:
      - Bins 1 & 8 (source/sink neighbours) always deviate strongly —
        they see the thermostat directly.
      - A systematic slope in the bulk bins signals non-uniform heat
        flux, which can happen at hetero-interfaces.
      - The interface bins (4|5) may show extra deviation if there is
        Kapitza resistance induced flux redistribution.
    """
    n = len(common)
    fig, axes = plt.subplots(2, n, figsize=(5 * n, 8),
                             gridspec_kw={'height_ratios': [1.4, 1.0]})
    if n == 1:
        axes = axes.reshape(2, 1)

    for col, f in enumerate(common):
        v    = results[f]
        bins = np.arange(1, N_BINS + 1)

        # Keep raw sign (negative flux direction) — matches your JP vs bin plot
        JP    = v['JP_mean']          # raw, negative values
        avg   = v['avg_jp']           # fit-bin mean, also negative
        dev   = JP - avg              # deviation from fit-bin mean

        # Block-averaged SEM — proper uncertainty on the time-averaged JP.
        # Split time series into N_BLOCKS non-overlapping blocks, average
        # each block, then take std of block means / sqrt(N_BLOCKS).
        # This accounts for autocorrelation; much smaller than per-step std.
        JP_ts    = v['JP_all']                              # (n_steps, N_BINS)
        N_BLOCKS = 10
        n_steps  = JP_ts.shape[0]
        blen     = n_steps // N_BLOCKS
        JP_trim  = JP_ts[:blen * N_BLOCKS, :]              # trim to exact multiple
        blk      = JP_trim.reshape(N_BLOCKS, blen, N_BINS)
        blk_mean = blk.mean(axis=1)                        # (N_BLOCKS, N_BINS)
        jp_std   = blk_mean.std(axis=0) / np.sqrt(N_BLOCKS)  # SEM of block means

        # ±1σ = spatial std of the fit-bin TIME-AVERAGED JP values
        # (how much the fit bins disagree with each other after averaging)
        # This is the meaningful "flatness" check — much smaller than temporal σ
        sigma_fit = JP[FIT_0IDX].std()

        # ── TOP: raw JP (same sign convention as plot_jp_vs_bin) ─
        ax_top = axes[0, col]
        ax_top.plot(bins, JP, 'o-', color='C0', lw=2, ms=8, label='All bins')
        ax_top.plot([bins[i] for i in FIT_0IDX], JP[FIT_0IDX],
                    'ro', ms=10, mfc='white', mew=2, label='Fit bins')
        ax_top.axhline(avg, color='C3', ls='--', lw=1.5,
                       label=f'avg={avg:.4f}')
        ax_top.axhline(0, color='k', lw=0.6, ls=':')
        ax_top.set_ylabel(f'{COMPONENT} (eV$^{{3/2}}$ amu$^{{-1/2}}$)', fontsize=12)
        ax_top.set_title(f'{f}  Q={v["Q_Wm2"]:.3e} W/m²', fontsize=11)
        ax_top.legend(fontsize=9)
        ax_top.grid(alpha=0.3)
        ax_top.set_xticks(bins)

        # shade source/sink neighbour bins (1 and 8) — always noisy
        for b in [1, N_BINS]:
            ax_top.axvspan(b - 0.45, b + 0.45, color='grey', alpha=0.10, zorder=0)

        # ── BOTTOM: deviation ─────────────────────────────────
        ax_bot = axes[1, col]

        # ±1σ band: spatial scatter of fit-bin means (true flux non-uniformity)
        ax_bot.axhspan(-sigma_fit, sigma_fit, color='C0', alpha=0.12,
                       label=f'±1σ fit bins ({sigma_fit:.5f})')

        # zero line
        ax_bot.axhline(0, color='k', lw=1.2, ls='-')

        # colour bars: above mean = red, below mean = blue
        colours = ['#c0392b' if d >= 0 else '#2471a3' for d in dev]
        ax_bot.bar(bins, dev, color=colours, width=0.55, alpha=0.75,
                   zorder=3, label='Deviation (bin − mean)')

        # error bars: temporal std per bin (uncertainty on the time average)
        ax_bot.errorbar(bins, dev, yerr=jp_std,
                        fmt='none', ecolor='k', elinewidth=1.2,
                        capsize=4, capthick=1.2, zorder=4,
                        label=f'Block SEM (n={N_BLOCKS})')

        # interface marker
        x_iface = 0.5 * (IFACE_LEFT_0IDX + 1 + IFACE_RIGHT_0IDX + 1)
        ax_bot.axvline(x_iface, color='purple', ls='--', lw=1.6,
                       alpha=0.8, label=f'Interface ({IFACE_LABEL})')

        # shade source/sink neighbour bins
        for b in [1, N_BINS]:
            ax_bot.axvspan(b - 0.45, b + 0.45, color='grey', alpha=0.10, zorder=0)
            ax_bot.text(b, dev[b - 1],
                        'src/snk\nneighbour',
                        ha='center',
                        va='bottom' if dev[b-1] >= 0 else 'top',
                        fontsize=7, color='grey')

        ax_bot.set_xlabel('Bin index', fontsize=13)
        ax_bot.set_ylabel('ΔJP  (bin − fit-mean)\n'
                          f'(eV$^{{3/2}}$ amu$^{{-1/2}}$)', fontsize=11)
        ax_bot.set_xticks(bins)
        ax_bot.grid(alpha=0.3)
        ax_bot.legend(fontsize=8, loc='upper right')

        # % deviation labels — skip bins 1 & 8 (always outliers)
        for b in bins:
            if b in [1, N_BINS]:
                continue
            pct  = dev[b-1] / abs(avg) * 100 if avg != 0 else 0
            yoff = dev[b-1] + np.sign(dev[b-1]) * abs(sigma_fit) * 0.25
            ax_bot.text(b, yoff, f'{pct:+.1f}%',
                        ha='center',
                        va='bottom' if dev[b-1] >= 0 else 'top',
                        fontsize=8, color='#333', fontweight='bold')

    plt.suptitle(f'JP vs Bin & Deviation from Mean  |  {COMPONENT}',
                 fontsize=13, y=1.01)
    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fname = f"JP_deviation_{COMPONENT}.{ext}"
        plt.savefig(fname, dpi=300, bbox_inches='tight')
        print(f"  Saved: {fname}")
    plt.show()
    plt.close()


# ============================================================
# PAPER-STYLE TEMPERATURE PROFILE PLOT
# ============================================================
def plot_paper_style(v):
    """Reproduce the paper-style temperature profile (mean dots + extended fit lines)."""
    f          = v['folder']
    x_all_nm   = v['x_all_nm']
    T_all      = v['T_all']
    T_mean     = v['T_mean']
    x_iface_nm = v['x_iface_m'] * 1e9
    L_nm       = v['L_nm']
    bin_w      = v['L_bin_nm']

    src_bounds = v.get('src_bounds_nm')
    snk_bounds = v.get('snk_bounds_nm')
    margin = 0.03 * L_nm

    if src_bounds and snk_bounds:
        x_left_edge        = -margin
        x_right_edge        = L_nm + margin
        source_band_left    = -margin
        source_band_right   = src_bounds[1]
        sink_band_left       = snk_bounds[0]
        sink_band_right      = L_nm + margin
    else:
        bin_w = v['L_bin_nm']
        x_left_edge  = max(0.0, x_all_nm[0]  - bin_w * 1.3)
        x_right_edge =          x_all_nm[-1] + bin_w * 1.3
        source_band_left  = x_left_edge
        source_band_right = x_all_nm[0]  + bin_w * 0.45
        sink_band_left    = x_all_nm[-1] - bin_w * 0.45
        sink_band_right   = x_right_edge

    C_SCATTER = '#f5a623'
    C_FIT     = '#1a3a8f'
    C_SOURCE  = '#c0392b'
    C_SINK    = '#2471a3'

    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    ax.axvspan(source_band_left, source_band_right, color=C_SOURCE, alpha=0.20, zorder=0)
    ax.axvspan(sink_band_left, sink_band_right,      color=C_SINK,   alpha=0.20, zorder=0)

    ax.plot(x_all_nm, T_mean, 'o', color=C_SCATTER,
            ms=9, mec='#b7590a', mew=1.2, zorder=3)

    xl = np.linspace(x_left_edge, x_iface_nm, 300)
    ax.plot(xl, v['grad_L'] * xl * 1e-9 + v['intc_L'],
            '-', color=C_FIT, lw=2.5, zorder=4)

    xr = np.linspace(x_iface_nm, x_right_edge, 300)
    ax.plot(xr, v['grad_R'] * xr * 1e-9 + v['intc_R'],
            '-', color=C_FIT, lw=2.5, zorder=4)

    ax.axvline(x_iface_nm, color=C_SINK, ls='--', lw=1.8, alpha=0.9, zorder=5)

    T_L = v['T_iface_L']
    T_R = v['T_iface_R']
    ax.annotate('', xy=(x_iface_nm, T_R), xytext=(x_iface_nm, T_L),
                arrowprops=dict(arrowstyle='<->', color='black', lw=2.0), zorder=6)
    ax.text(x_iface_nm + 0.4, 0.5 * (T_L + T_R),
            f'$\\Delta T$ = {v["delta_T_iface"]:.2f} K',
            color='black', fontsize=16, va='center', fontweight='bold',
            fontfamily='serif')

    y_mid = 0.5 * (T_mean.min() + T_mean.max())
    x_source_label = 0.5 * (source_band_left + source_band_right)
    x_sink_label   = 0.5 * (sink_band_left + sink_band_right)
    ax.text(x_source_label, y_mid, 'Heat\nsource',
            color=C_SOURCE, fontsize=13, fontweight='bold',
            va='center', ha='center', rotation=90, fontfamily='serif',
            zorder=7, clip_on=False)
    ax.text(x_sink_label, y_mid, 'Heat\nsink',
            color=C_SINK, fontsize=13, fontweight='bold',
            va='center', ha='center', rotation=90, fontfamily='serif',
            zorder=7, clip_on=False)

    ax.text(0.25, 0.95, 'Bi$_2$Te$_3$', transform=ax.transAxes,
            ha='center', va='top', color=C_SOURCE, fontsize=18,
            fontweight='bold', fontfamily='serif')
    ax.text(0.75, 0.95, 'Sb$_2$Te$_3$', transform=ax.transAxes,
            ha='center', va='top', color=C_SINK, fontsize=18,
            fontweight='bold', fontfamily='serif')

    domain_w      = x_right_edge - x_left_edge
    x_arrow_end   = x_iface_nm - 0.10 * domain_w
    x_arrow_start = x_iface_nm - 0.32 * domain_w
    x_text        = x_iface_nm - 0.21 * domain_w
    arrow_y = 0.5 * (v['T_iface_L'] + v['T_iface_R'])
    ax.annotate('', xy=(x_arrow_end, arrow_y),
                xytext=(x_arrow_start, arrow_y),
                arrowprops=dict(arrowstyle='->', color='#333', lw=2.0), zorder=6)
    ax.text(x_text, arrow_y + (T_mean.max() - T_mean.min()) * 0.04,
            'Heat flux', fontsize=14, color='#333', ha='center', fontfamily='serif',
            zorder=6)

    ax.set_xlim(x_left_edge - 0.15 * bin_w, x_right_edge + 0.15 * bin_w)
    ax.set_xlabel('Position [nm]')
    ax.set_ylabel('Temperature [K]')
    ax.set_title('')
    ax.grid(alpha=0.25)
    apply_style(ax)
    plt.tight_layout()

    for ext in ('pdf', 'png'):
        fname = f"T_profile_paper_{f}.{ext}"
        plt.savefig(fname, dpi=300, bbox_inches='tight')
        print(f"  Saved: {fname}")
    plt.close()


# ============================================================
# MAIN
# ============================================================
if os.path.exists('compute.out'):
    FOLDERS     = ['.']
    SINGLE_MODE = True
else:
    FOLDERS = sorted([d for d in os.listdir('.')
                      if d.startswith(('L_', 'SI_')) and os.path.isdir(d)])
    SINGLE_MODE = False
    if not FOLDERS:
        print("No compute.out found here and no L_xxxx folders found!"); sys.exit(1)

print(f"Mode      : {'Single folder (cwd)' if SINGLE_MODE else 'Multi-folder'}")
print(f"Folders   : {FOLDERS}")
print(f"N_BINS={N_BINS}  REMOVE_N={REMOVE_N}  FIT_GIDS={FIT_GIDS}")
print(f"COMPONENT : {COMPONENT}")
print(f"Interface : {IFACE_LABEL}  (bins {IFACE_LEFT_0IDX+1}|{IFACE_RIGHT_0IDX+1})")

results = {}
for folder in FOLDERS:
    label = os.path.basename(os.path.abspath(folder)) if folder == '.' else folder
    print(f"\n{'='*60}\n  {label}\n{'='*60}")
    r = compute_folder(folder)
    if r: results[label] = r

if not results:
    print("\nNo valid folders!"); sys.exit(0)

common = sorted(results.keys())
common = sorted(common, key=lambda f: results[f]['L_nm'])

L_arr  = np.array([results[f]['L_nm']  for f in common])
k_arr  = np.array([results[f]['kappa'] for f in common])
RK_arr = np.array([results[f]['R_K']   for f in common])
GK_arr = np.array([results[f]['G_K']   for f in common])
dT_arr = np.array([results[f]['delta_T_iface'] for f in common])

# ============================================================
# TABLE
# ============================================================
W = 80
print(f"\n{'='*W}")
print(f"  INTERFACE THERMAL RESISTANCE  |  {IFACE_LABEL}  |  {COMPONENT}")
print(f"{'='*W}")
print(f"  {'Folder':<30} {'L(nm)':>7} {'kappa':>9} {'dT_IF(K)':>10} "
      f"{'R_K(m2K/W)':>14} {'R_K(1e-9)':>12} {'G_K(MW/m2K)':>13}")
print(f"  {'-'*78}")
for f in common:
    v = results[f]
    print(f"  {f:<30} {v['L_nm']:>7.1f} {v['kappa']:>9.4f} "
          f"{v['delta_T_iface']:>10.4f} {v['R_K']:>14.4e} "
          f"{v['R_K_nm2KW']:>13.4f} {v['G_K']/1e6:>13.4f}")
if len(results) >= 2:
    print(f"\n  Mean R_K = {RK_arr.mean():.4e} m2K/W  ({RK_arr.mean()*1e9:.4f} x1e-9 m2K/W)")
    print(f"  Std  R_K = {RK_arr.std():.4e} m2K/W")
    print(f"  Mean G_K = {GK_arr.mean()/1e6:.4f} MW/m2K")

# FSC
if len(results) >= 2:
    iL = 1.0 / L_arr; ik = 1.0 / k_arr
    cf = np.polyfit(iL, ik, 1)
    kb = 1.0 / cf[1]
    r2 = 1 - np.sum((ik - np.polyval(cf, iL))**2) / np.sum((ik - ik.mean())**2)
    print(f"\n  FSC: kappa_inf = {kb:.4f} W/mK  (R2={r2:.4f})")
else:
    kb, r2, cf, iL, ik = None, None, None, None, None
print(f"{'='*W}")

# ============================================================
# PLOT 1a: JP vs Bin
# ============================================================
plot_jp_vs_bin(common, results)

# ============================================================
# PLOT 1b: JP Deviation from Mean  ← NEW
# ============================================================
plot_jp_deviation(common, results)

# ============================================================
# PLOT 2: Paper-style temperature profiles (one per folder)
# ============================================================
for f in common:
    plot_paper_style(results[f])

# ============================================================
# PLOT 3: Summary (multi-folder) — R_K, dT_IF, FSC
# ============================================================
def savefig(name):
    plt.savefig(f"{name}.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{name}.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {name}.pdf / .png")

if len(results) >= 2:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    ax = axes[0]
    ax.plot(L_arr, RK_arr * 1e9, 's-', color='purple', lw=2, ms=10, mfc='white', mew=2)
    ax.axhline(RK_arr.mean() * 1e9, color='purple', ls='--', lw=1.5,
               label=f'mean={RK_arr.mean()*1e9:.3f} x1e-9 m2K/W')
    ax.set_xlabel('L (nm)', fontsize=13); ax.set_ylabel('R_K (x1e-9 m2K/W)', fontsize=13)
    ax.set_title('Kapitza Resistance vs L', fontsize=11)
    ax.legend(fontsize=10); ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(L_arr, dT_arr, '^-', color='darkorange', lw=2, ms=10, mfc='white', mew=2)
    ax.set_xlabel('L (nm)', fontsize=13); ax.set_ylabel('dT_interface (K)', fontsize=13)
    ax.set_title('Interface Temperature Jump vs L', fontsize=11); ax.grid(alpha=0.3)

    ax = axes[2]
    x_ext2 = np.linspace(0, iL.max() * 1.2, 300)
    ax.plot(iL, ik, 'o', color='C0', ms=11, mfc='white', mew=2, zorder=5, label='NEMD')
    ax.plot(x_ext2, np.polyval(cf, x_ext2), 'k--', lw=2, label=f'Fit R2={r2:.4f}')
    ax.plot(0, cf[1], 'r*', ms=15, zorder=6, label=f'kappa_inf={kb:.3f} W/mK')
    ax.set_xlim(left=-0.002)
    ax.set_xlabel('1/L (nm-1)', fontsize=13); ax.set_ylabel('1/kappa (mK/W)', fontsize=13)
    ax.set_title('Finite-size correction', fontsize=11); ax.legend(fontsize=10)
    plt.tight_layout()
    savefig(f"RK_summary_{COMPONENT}")

print(f"\n{'='*W}")
print(f"  FINAL  |  {IFACE_LABEL}  |  {COMPONENT}")
print(f"{'='*W}")
for f in common:
    v = results[f]
    print(f"  {f}: R_K={v['R_K']:.4e} m2K/W  ({v['R_K_nm2KW']:.4f} x1e-9 m2K/W)  G_K={v['G_K']/1e6:.4f} MW/m2K")
if len(results) >= 2:
    print(f"\n  Mean R_K = {RK_arr.mean():.4e} +/- {RK_arr.std():.4e} m2K/W")
    print(f"  Mean G_K = {GK_arr.mean()/1e6:.4f} +/- {GK_arr.std()/1e6:.4f} MW/m2K")
    print(f"  kappa_inf = {kb:.4f} W/mK  (R2={r2:.4f})")
print(f"{'='*W}")
