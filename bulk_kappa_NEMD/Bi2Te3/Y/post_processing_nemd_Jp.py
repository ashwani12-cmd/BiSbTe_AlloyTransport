import numpy as np
import matplotlib.pyplot as plt
import os, sys, re

# ============================================================
# USER SETTINGS
# ============================================================
N_BINS    = 8
REMOVE_N  = 2
COMPONENT = 'JPy'   # JPx | JPy | JPz

# None = ALL rows | e.g. 2500 = last 2500 rows only
JP_STEPS  = None
# ============================================================

dt_fs           = 1.0
sample_interval = 10
output_interval = 100
FRAC_WALL       = 0.02
FRAC_SRC        = 0.18
FRAC_SNK        = 0.18

time_per_row_ns = dt_fs * sample_interval * output_interval * 1e-6
eV_to_J         = 1.602176634e-19
amu_to_kg       = 1.66053906660e-27
JP_conv         = (eV_to_J ** 1.5) * (amu_to_kg ** -0.5)

FIT_0IDX = list(range(REMOVE_N, N_BINS - REMOVE_N))   # [2,3,4,5]
FIT_GIDS = [i + 1 for i in FIT_0IDX]                  # [3,4,5,6]

# ============================================================
# HELPERS
# ============================================================

def parse_setup(folder):
    L, A = None, None
    try:
        with open(os.path.join(folder, 'nemd_setup.txt')) as f:
            for line in f:
                s = line.strip()
                if 'L_transport' in s and ':' in s and 'A' in s:
                    try: L = float(s.split(':')[1].split('A')[0].strip())
                    except: pass
                if 'A_cross' in s and ':' in s and 'A2' in s:
                    try: A = float(s.split(':')[1].split('A2')[0].strip())
                    except: pass
    except: pass
    return L, A


def read_bin_centers_nm(folder):
    centers = {}
    try:
        with open(os.path.join(folder, 'nemd_setup.txt')) as f:
            in_sec = False
            for line in f:
                s = line.strip()
                if 'BIN CENTERS' in s:
                    in_sec = True; continue
                if in_sec:
                    if s.startswith('===') or s == '': break
                    nums = re.findall(r"[-+]?\d*\.\d+|\d+", s)
                    if len(nums) >= 4:
                        gid = int(nums[0])
                        centers[gid] = float(nums[-2])   # nm
    except: pass
    return centers


def fit_gradient(x_m, T_vals):
    coeffs = np.polyfit(x_m, T_vals, 1)
    T_pred = np.polyval(coeffs, x_m)
    ss_res = np.sum((T_vals - T_pred)**2)
    ss_tot = np.sum((T_vals - T_vals.mean())**2)
    R2 = 1 - ss_res/ss_tot if ss_tot > 0 else 1.0
    return coeffs[0], coeffs[1], R2


def get_jp_col_start(ncols):
    comp      = COMPONENT.lower()
    block_idx = {'jpx':0,'jpy':1,'jpz':2,'jkx':3,'jky':4,'jkz':5}[comp]
    T_end     = 1 + N_BINS
    jp_start  = T_end + block_idx*(N_BINS+1) + 1
    return jp_start


# ============================================================
# COMPUTE ONE FOLDER
# ============================================================

def compute_folder(folder):
    fpath = os.path.join(folder, 'compute.out')
    if not os.path.exists(fpath):
        print(f"  SKIP {folder} — compute.out not found")
        return None

    L_ang, A_ang2 = parse_setup(folder)
    if L_ang is None:
        print(f"  SKIP {folder} — nemd_setup.txt parse failed")
        return None

    data         = np.loadtxt(fpath)
    ndata, ncols = data.shape
    t_ns         = np.arange(ndata) * time_per_row_ns

    # ── Step range ───────────────────────────────────────────
    jp_start_row = 0 if JP_STEPS is None else max(0, ndata - JP_STEPS)
    steps_jp     = ndata - jp_start_row

    A_m2      = A_ang2 * 1e-20
    L_nm      = L_ang  * 0.1
    A_nm2     = A_ang2 * 1e-2
    L_bin_ang = L_ang * (1.0 - 2*FRAC_WALL - FRAC_SRC - FRAC_SNK) / N_BINS
    L_bin_m   = L_bin_ang * 1e-10
    L_bin_nm  = L_bin_ang * 0.1
    Vol_nm3   = A_nm2 * L_bin_nm

    # ── Bin centers ───────────────────────────────────────────
    centers_nm = read_bin_centers_nm(folder)
    x_fit_m    = np.array([centers_nm[g] for g in FIT_GIDS]) * 1e-9

    # ── Temperature (jp_start_row..end) ─────────────────────
    T_all  = data[jp_start_row:, 1:N_BINS+1]
    T_mean = np.mean(T_all, axis=0)
    T_fit  = T_mean[FIT_0IDX]
    grad_T, intercept, R2_T = fit_gradient(x_fit_m, T_fit)
    deltaT = T_mean[0] - T_mean[N_BINS-2]

    # ── JP (jp_start_row..end) ───────────────────────────────
    jp_col   = get_jp_col_start(ncols)
    JP_all   = data[jp_start_row:, jp_col:jp_col+N_BINS]
    JP_mean  = np.mean(JP_all, axis=0)
    avg_jp   = np.mean(JP_mean[FIT_0IDX])
    Q_Wm2    = (avg_jp * JP_conv) / (A_m2 * L_bin_m)
    kappa    = -Q_Wm2 / grad_T

    print(f"\n  {folder}  |  L={L_nm:.2f}nm  A={A_nm2:.3f}nm²  "
          f"L_bin={L_bin_nm:.4f}nm  Vol={Vol_nm3:.3f}nm³")
    print(f"  Steps: {steps_jp} (rows {jp_start_row}–{ndata})  "
          f"{'ALL' if JP_STEPS is None else f'last {JP_STEPS}'}")
    print(f"  avg_JP={avg_jp:.6f}  Q={Q_Wm2:.4e} W/m²")
    print(f"  dT/dx={grad_T:.4e} K/m  R²={R2_T:.6f}  ΔT={deltaT:.3f}K")
    print(f"  κ ({COMPONENT}) = {kappa:.4f} W/mK")

    return dict(
        folder=folder, L_nm=L_nm, A_nm2=A_nm2, A_m2=A_m2,
        L_bin_nm=L_bin_nm, Vol_nm3=Vol_nm3,
        ndata=ndata, steps_jp=steps_jp,
        T_mean=T_mean, T_fit=T_fit,
        x_fit_m=x_fit_m, centers_nm=centers_nm,
        grad_T=grad_T, intercept=intercept, R2_T=R2_T, deltaT=deltaT,
        JP_mean=JP_mean, JP_all=JP_all, avg_jp=avg_jp,
        Q_Wm2=Q_Wm2, kappa=kappa, t_ns=t_ns,
    )


# ============================================================
# MAIN
# ============================================================
FOLDERS = sorted([d for d in os.listdir('.')
                  if d.startswith('L_') and os.path.isdir(d)])
if not FOLDERS:
    print("No L_xxxx folders found!"); sys.exit(1)

print(f"Folders   : {FOLDERS}")
print(f"N_BINS={N_BINS}  REMOVE_N={REMOVE_N}  "
      f"FIT_0IDX={FIT_0IDX}  FIT_GIDS={FIT_GIDS}")
print(f"COMPONENT : {COMPONENT}")
print(f"JP_STEPS  : {'ALL' if JP_STEPS is None else JP_STEPS}")

results = {}
for folder in FOLDERS:
    print(f"\n{'='*60}\n  {folder}\n{'='*60}")
    r = compute_folder(folder)
    if r: results[folder] = r

if len(results) < 2:
    print("\nNeed ≥2 folders!"); sys.exit(0)

# ============================================================
# SORT
# ============================================================
common = sorted(results.keys())
L_arr  = np.array([results[f]['L_nm']  for f in common])
k_arr  = np.array([results[f]['kappa'] for f in common])

# ============================================================
# FINITE-SIZE CORRECTION
# ============================================================
iL  = 1.0 / L_arr
ik  = 1.0 / k_arr
cf  = np.polyfit(iL, ik, 1)
kb  = 1.0 / cf[1]
r2  = 1 - np.sum((ik-np.polyval(cf,iL))**2)/np.sum((ik-ik.mean())**2)

# ============================================================
# TABLE 1 — Per-length details
# ============================================================
W = 75
print(f"\n{'='*W}")
print(f"  TABLE 1: Per-length details  |  {COMPONENT} method")
print(f"  Fit bins (0-idx): {FIT_0IDX} = groups {FIT_GIDS}")
print(f"  Steps: {'ALL' if JP_STEPS is None else JP_STEPS}")
print(f"{'='*W}")

# Section A: Geometry + gradient
print(f"\n  -- A: Geometry & Temperature Gradient --")
print(f"  {'Folder':<10} {'L(nm)':>6} {'A(nm²)':>9} {'L_bin(nm)':>10} "
      f"{'Vol(nm³)':>10} {'ΔT(K)':>7} {'dT/dx(K/m)':>13} {'R²':>8}")
print(f"  {'─'*77}")
for f in common:
    v = results[f]
    print(f"  {f:<10} {v['L_nm']:>6.1f} {v['A_nm2']:>9.3f} "
          f"{v['L_bin_nm']:>10.4f} {v['Vol_nm3']:>10.3f} "
          f"{v['deltaT']:>7.3f} {v['grad_T']:>13.4e} {v['R2_T']:>8.6f}")

# Section B: Heat flux + kappa
print(f"\n  -- B: Heat Flux & Thermal Conductivity --")
print(f"  {'Folder':<10} {'avg_JP':>12} {'Q(W/m²)':>13} "
      f"{'κ(W/mK)':>9} {'Steps':>7}")
print(f"  {'─'*55}")
for f in common:
    v = results[f]
    print(f"  {f:<10} {v['avg_jp']:>12.6f} {v['Q_Wm2']:>13.4e} "
          f"{v['kappa']:>9.4f} {v['steps_jp']:>7}")
print(f"  {'─'*55}")
print(f"  {'∞':<10} {'—':>12} {'—':>13} {kb:>9.4f}")
print(f"{'='*W}")

# ============================================================
# TABLE 2 — FSC fitting
# ============================================================
# Detect missing standard lengths
std_A     = [250, 500, 750, 1000]
present_A = [int(round(l*10)) for l in L_arr]
missing_A = [x for x in std_A if x not in present_A]

print(f"\n{'='*W}")
print(f"  TABLE 2: Finite-Size Correction  (1/κ vs 1/L)  |  {COMPONENT}")
print(f"  Lengths used : {[f'{l:.1f}nm' for l in L_arr]}")
if missing_A:
    print(f"  ⚠️  Missing  : {missing_A} Å — skipped (fit uses available only)")
print(f"{'='*W}")

print(f"\n  {'L(nm)':>8} {'1/L(nm⁻¹)':>12} {'κ(W/mK)':>10} "
      f"{'1/κ':>10} {'1/κ_fit':>10}")
print(f"  {'─'*55}")
ik_fit = np.polyval(cf, iL)
for i in range(len(L_arr)):
    print(f"  {L_arr[i]:>8.2f} {iL[i]:>12.6f} {k_arr[i]:>10.4f} "
          f"{ik[i]:>10.6f} {ik_fit[i]:>10.6f}")
print(f"  {'─'*55}")
print(f"\n  Fit:      1/κ = {cf[0]:.6f} × (1/L) + {cf[1]:.6f}")
print(f"  κ∞      = 1 / {cf[1]:.6f} = {kb:.4f} W/mK")
print(f"  R²      = {r2:.6f}  {'✅' if r2>0.95 else '⚠️ R² low'}")
print(f"  Slope   = {cf[0]:.4f} nm·mK/W")
print(f"  MFP est = slope × κ∞ ≈ {cf[0]*kb:.2f} nm  (rough)")
print(f"  Points  = {len(L_arr)}  (250→{int(L_arr[-1]*10)} Å)")
print(f"\n  🔥 κ∞ ({COMPONENT}) = {kb:.4f} W/mK  (R²={r2:.4f})")
print(f"{'='*W}")


# ============================================================
# PLOTS
# ============================================================
def savefig(name):
    plt.savefig(f"{name}.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {name}.pdf")

n = len(common)

# ── Plot 1: T profiles ───────────────────────────────────────
fig, axes = plt.subplots(1, n, figsize=(5*n, 4.5))
if n == 1: axes = [axes]
for ax, f in zip(axes, common):
    v        = results[f]
    x_all_nm = np.array([v['centers_nm'][g] for g in range(1, N_BINS+1)])
    x_fit_nm = v['x_fit_m'] * 1e9
    x_line   = np.linspace(x_fit_nm[0], x_fit_nm[-1], 100)
    T_line   = v['grad_T'] * x_line * 1e-9 + v['intercept']
    ax.plot(x_all_nm, v['T_mean'], 'o-', color='gray', alpha=0.5, ms=6, label='All bins')
    ax.plot(x_fit_nm, v['T_fit'],  'ro', ms=9, mfc='white', mew=2, label='Fit bins')
    ax.plot(x_line, T_line, 'k--', lw=2, label='Linear fit')
    ax.set_xlabel('Position (nm)', fontsize=13)
    ax.set_ylabel('$T$ (K)', fontsize=13)
    ax.set_title(f'{f}\nκ={v["kappa"]:.3f} W/mK  dT/dx={v["grad_T"]:.3e}  '
                 f'R²={v["R2_T"]:.4f}', fontsize=10)
    ax.legend(fontsize=8)
plt.tight_layout()
savefig(f"temperature_profiles_{COMPONENT}")

'''
# ── Plot 2: JP vs time ───────────────────────────────────────
fig, axes = plt.subplots(1, n, figsize=(7*n, 4.5))
if n == 1: axes = [axes]
colors = plt.cm.Blues(np.linspace(0.3, 0.9, N_BINS))
for ax, f in zip(axes, common):
    v = results[f]
    for i in range(N_BINS):
        lw   = 2.0 if i in FIT_0IDX else 0.7
        alph = 1.0 if i in FIT_0IDX else 0.3
        ax.plot(v['t_ns'], v['JP_all'][:, i],
                lw=lw, alpha=alph, color=colors[i],
                label=f'G{i+1}{"(fit)" if i in FIT_0IDX else ""}')
    ax.axhline(v['avg_jp'], color='C3', ls='--', lw=1.5,
               label=f'avg={v["avg_jp"]:.4f}')
    ax.axhline(0, color='k', lw=0.6, ls=':')
    ax.set_xlabel('$t$ (ns)', fontsize=13)
    ax.set_ylabel(f'{COMPONENT} (eV$^{{3/2}}$ amu$^{{-1/2}}$)', fontsize=12)
    ax.set_title(f'{f}  avg(fit)={v["avg_jp"]:.4f}  '
                 f'Q={v["Q_Wm2"]:.3e} W/m²', fontsize=10)
    ax.legend(fontsize=7, ncol=2)
plt.tight_layout()
savefig(f"JP_vs_time_{COMPONENT}")
'''

# ── Plot 3: JP vs bin ────────────────────────────────────────
fig, axes = plt.subplots(1, n, figsize=(5*n, 4))
if n == 1: axes = [axes]
for ax, f in zip(axes, common):
    v    = results[f]
    bins = np.arange(1, N_BINS+1)
    ax.plot(bins, -v['JP_mean'], 'o-', color='C0', lw=2, ms=8, label='All bins')
    ax.plot([bins[i] for i in FIT_0IDX], -v['JP_mean'][FIT_0IDX],
            'ro', ms=10, mfc='white', mew=2, label='Fit bins')
    ax.axhline(-v['avg_jp'], color='C3', ls='--', lw=1.5,
               label=f'avg={-v["avg_jp"]:.4f}')
    ax.axhline(0, color='k', lw=0.6, ls=':')
    ax.set_xlabel('Bin index', fontsize=13)
    ax.set_ylabel(f'{COMPONENT} (eV$^{{3/2}}$ amu$^{{-1/2}}$)', fontsize=12)
    ax.set_title(f'{f}  Q={-v["Q_Wm2"]:.3e} W/m²', fontsize=11)
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
plt.show()
savefig(f"JP_vs_bin_{COMPONENT}")

# ── Plot 4: κ vs L ───────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
ax = axes[0]
ax.plot(L_arr, k_arr, 'o-', color='C0', lw=2, ms=10, mfc='white', mew=2)
ax.axhline(kb, color='red', ls='--', lw=1.5, label=f'κ∞={kb:.3f} W/mK')
ax.set_xlabel('$L$ (nm)', fontsize=14)
ax.set_ylabel(r'$\kappa$ (W m$^{-1}$ K$^{-1}$)', fontsize=14)
ax.set_title(f'Size-dependent κ  ({COMPONENT})', fontsize=13)
ax.legend(fontsize=12); ax.grid(alpha=0.4)

ax = axes[1]
x_ext = np.linspace(0, iL.max()*1.2, 300)
ax.plot(iL, ik, 'o', color='C0', ms=11, mfc='white', mew=2, zorder=5, label='NEMD')
ax.plot(x_ext, np.polyval(cf, x_ext), 'k--', lw=2,
        label=f'Linear fit  R²={r2:.4f}')
ax.plot(0, cf[1], 'r*', ms=15, zorder=6, label=f'κ∞={kb:.3f} W/mK')
ax.set_xlim(left=-0.002)
ax.set_xlabel('$1/L$ (nm$^{-1}$)', fontsize=14)
ax.set_ylabel(r'$1/\kappa$ (m K W$^{-1}$)', fontsize=14)
ax.set_title(f'Finite-size correction  ({COMPONENT})', fontsize=13)
ax.legend(fontsize=11)
plt.tight_layout()
plt.show()
savefig(f"finite_size_{COMPONENT}")

print(f"\n{'='*W}")
print(f"  🔥 κ∞ ({COMPONENT}) = {kb:.4f} W/mK  (R²={r2:.4f})")
print(f"{'='*W}")
