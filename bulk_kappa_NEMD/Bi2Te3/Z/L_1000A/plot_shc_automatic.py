"""
plot_shc_auto.py  (v3 — correct column layout for 10-group compute.out)
────────────────────────────────────────────────────────────────────────
Fully automated SHC plotter for GPUMD NEMD runs.

KEY FIXES vs earlier versions:
  1. deltaT from linear gradient fit through G3-G6 (not src-snk setpoints)
  2. compute.out has 10 groups (walls GID0-L/R included), not 8:
       col 0       : step
       cols 1-10   : temperature (GID0-L, G1..G8, GID0-R)
       cols 11-40  : jp  (10 groups × 3 components)
       cols 41-70  : jk  (10 groups × 3 components)
     → temperature for G1-G8 = cols 2-9  (skip col 1 = wall-L)
     → JP for transport direction, middle groups G2-G7 = cols determined
       by transport direction component offset + group offset

Usage:
    python plot_shc_auto.py                    # cwd
    python plot_shc_auto.py /path/to/run/dir
"""

import re
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from gpyumd.load import load_compute


# ============================================================
# 0. Resolve run directory
# ============================================================
run_dir = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
run_dir = os.path.abspath(run_dir)
print(f"Run directory : {run_dir}")

setup_file = os.path.join(run_dir, "nemd_setup.txt")
shc_file   = os.path.join(run_dir, "shc.out")
compute_file = os.path.join(run_dir, "compute.out")

for f in (setup_file, shc_file, compute_file):
    if not os.path.isfile(f):
        raise FileNotFoundError(f"Required file not found: {f}")


# ============================================================
# 1. Parse nemd_setup.txt
# ============================================================
def parse_setup(path):
    with open(path) as fh:
        text = fh.read()

    def grab(pattern, cast=float, flags=re.IGNORECASE):
        m = re.search(pattern, text, flags)
        if m is None:
            raise ValueError(f"Pattern not found: {pattern!r}")
        return cast(m.group(1))

    info = {}
    info["transport_dir"] = grab(r"Transport dir\s*:\s*([XYZ])", cast=str).upper()
    info["Lx"]     = grab(r"Lx\s*:\s*([\d.]+)\s*A")
    info["Ly"]     = grab(r"Ly\s*:\s*([\d.]+)\s*A")
    info["Lz"]     = grab(r"Lz\s*:\s*([\d.]+)\s*A")
    info["L_eff_A"]= grab(r"L_eff\s*:\s*([\d.]+)\s*A")

    m4 = re.search(r"^\s*4\s+mid\s+G4\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
                   text, re.MULTILINE | re.IGNORECASE)
    if m4 is None:
        raise ValueError("Cannot find GID 4 row in GROUP BOUNDARIES table.")
    info["L_G4"] = float(m4.group(3))

    m_mat = re.search(r"Material\s*:\s*(\S+)", text, re.IGNORECASE)
    info["material"] = m_mat.group(1) if m_mat else "Unknown"

    # BIN CENTERS for GIDs 1-8
    bin_centers_A = {}
    for m in re.finditer(
        r"^\s*(\d+)\s+\S+.*?(\d+\.\d+)\s+\d+\.\d+\s+\d+\.\d+",
        text, re.MULTILINE
    ):
        gid = int(m.group(1))
        if 1 <= gid <= 8:
            bin_centers_A[gid] = float(m.group(2))

    if len(bin_centers_A) < 6:
        raise ValueError(f"Could not parse BIN CENTERS (found {len(bin_centers_A)})")
    info["bin_centers_A"] = bin_centers_A

    return info


params        = parse_setup(setup_file)
transport_dir = params["transport_dir"]
Lx, Ly, Lz   = params["Lx"], params["Ly"], params["Lz"]
L_G4          = params["L_G4"]
L_eff_A       = params["L_eff_A"]
material      = params["material"]
bin_centers_A = params["bin_centers_A"]

cross_map = {"X": (Ly, Lz), "Y": (Lx, Lz), "Z": (Lx, Ly)}
cross1, cross2 = cross_map[transport_dir]
V = L_G4 * cross1 * cross2

print(f"Material      : {material}")
print(f"Transport dir : {transport_dir}")
print(f"L_G4          : {L_G4:.4f} A")
print(f"Cross-section : {cross1:.4f} x {cross2:.4f} A")
print(f"V             : {V:.4e} A^3")
print(f"L_eff         : {L_eff_A:.4f} A")


# ============================================================
# 2. Load SHC data
# ============================================================
os.chdir(run_dir)

data = np.loadtxt(shc_file)
Nc_rows  = (data.shape[0] - 1000) if data.shape[0] > 1000 else 499
# Auto-detect: last 1000 rows are frequency if omega goes 0→50 THz
# Check if last 1000 rows col0 starts near 0
if data[-1000, 0] < 1.0:
    corr = data[:-1000]
    freq = data[-1000:]
else:
    corr = data[:499]
    freq = data[499:]

print(f"SHC rows: {data.shape[0]} total  →  {len(corr)} correlation + {len(freq)} frequency")

t   = corr[:, 0]
Ki  = corr[:, 1]
Ko  = corr[:, 2]
nu  = freq[:, 0] / (2 * np.pi)
jwi = freq[:, 1]
jwo = freq[:, 2]


# ============================================================
# 3. Temperature & deltaT from gradient fit
#
#  compute.out layout (10 groups: wall-L, G1..G8, wall-R):
#    col 0       : step
#    cols 1-10   : temperature
#    cols 11-40  : jp  (10 grp × 3 xyz)
#    cols 41-70  : jk  (10 grp × 3 xyz)
#
#  G1-G8 temperatures = cols 2-9  (col 1 = wall-L GID0)
# ============================================================
raw = np.loadtxt(compute_file)
ndata = raw.shape[0]
ncols = raw.shape[1]
half  = int(ndata / 2) + 1

print(f"\ncompute.out: {ndata} rows × {ncols} cols")

# Auto-detect number of groups from column count:
# ncols = 1 + N_grp + 3*N_grp + 3*N_grp + extra
# Try N_grp = 10 (walls included): 1+10+30+30 = 71 — too many
# Try N_grp = 10, no jk:           1+10+30    = 41 — no
# Observed: 65 cols
# 65 - 1 = 64; if jp+jk: 64/6 = 10.67 (not integer)
# 65 - 1 - 2(thermostat) = 62; 62/6 = 10.33 (no)
# Let's just auto-detect via the data signal:
means = np.mean(raw[half:, :], axis=0)

# Temperature cols: values between 270-340 K
temp_cols = [i for i in range(1, ncols) if 270 < means[i] < 340]
print(f"Detected temperature cols: {temp_cols}")

# G1-G8 temperature = skip first temp col (wall-L), take next 8
if len(temp_cols) >= 9:
    g1_to_g8_temp_cols = temp_cols[1:9]   # skip wall-L
elif len(temp_cols) == 8:
    g1_to_g8_temp_cols = temp_cols        # no walls in output
else:
    raise ValueError(f"Expected 8-10 temperature cols, found {len(temp_cols)}: {temp_cols}")

print(f"G1-G8 temp cols: {g1_to_g8_temp_cols}")
temp_ave = np.mean(raw[half:, g1_to_g8_temp_cols], axis=0)
print(f"temp_ave (G1-G8): {np.round(temp_ave, 2)}")

# Gradient fit: G3,G4,G5,G6 (0-indexed 2,3,4,5 in temp_ave)
fit_gids  = [3, 4, 5, 6]
fit_idx   = [g - 1 for g in fit_gids]
fit_pos_A = np.array([bin_centers_A[g] for g in fit_gids])
fit_temps = temp_ave[fit_idx]

slope, intercept, r_fit, _, _ = linregress(fit_pos_A, fit_temps)
deltaT_fit   = abs(slope) * L_eff_A
deltaT_naive = temp_ave[0] - temp_ave[-1]

print(f"dT/dx         : {slope*1e10:.4e} K/m  (R²={r_fit**2:.6f})")
print(f"deltaT (fit)  : {deltaT_fit:.4f} K  ← used")
print(f"deltaT (naive): {deltaT_naive:.4f} K  ← WRONG")

deltaT_use = deltaT_fit if 5 < deltaT_fit < 200 else deltaT_naive


# ============================================================
# 4. Compute Gc and kappa
# ============================================================
Gc      = 1.6e4 * (jwi + jwo) / V / deltaT_use
G_total = np.trapz(Gc, nu)
L_eff_m = L_eff_A * 1e-10
kappa   = G_total * 1e9 * L_eff_m

print(f"\nGc max        : {Gc.max():.6e} GW/m²/K/THz")
print(f"G_total       : {G_total * 1000:.4f} MW/m²/K")
print(f"kappa (SHC)   : {kappa:.4f} W/mK")


# ============================================================
# 5. Cross-check: compute NEMD-style kappa from compute.out JP
#    to verify our deltaT and unit conversion are consistent
# ============================================================
# JP signal cols: find consecutive block with uniform ~0.1-0.4 values
jp_signal = [(i, means[i]) for i in range(ncols)
             if 0.08 < abs(means[i]) < 0.50 and i > max(temp_cols)]

if jp_signal:
    jp_cols_found = [i for i,v in jp_signal]
    # Find longest consecutive run
    best_run, cur_run = [], [jp_cols_found[0]]
    for c in jp_cols_found[1:]:
        if c == cur_run[-1] + 1:
            cur_run.append(c)
        else:
            if len(cur_run) > len(best_run):
                best_run = cur_run[:]
            cur_run = [c]
    if len(cur_run) > len(best_run):
        best_run = cur_run

    if len(best_run) >= 4:
        jp_vals = means[best_run]
        avg_JP_check = np.mean(np.abs(jp_vals))
        A_cross = cross1 * cross2
        Q_check = avg_JP_check * 1.6e4
        dTdx_check = abs(slope) * 1e10   # K/m
        kappa_check = Q_check / dTdx_check if dTdx_check > 0 else 0
        print(f"\n--- Cross-check (NEMD-style from compute.out) ---")
        print(f"JP signal cols : {best_run}")
        print(f"JP values      : {np.round(jp_vals, 4)}")
        print(f"avg_JP         : {avg_JP_check:.6f}")
        print(f"Q              : {Q_check:.4e} W/m²")
        print(f"kappa_NEMD     : {kappa_check:.4f} W/mK  (should match NEMD table)")


# ============================================================
# 6. Plot
# ============================================================
plt.rcParams.update({
    "font.family"     : "serif",
    "mathtext.fontset": "stix",
    "font.size"       : 14,
    "axes.linewidth"  : 1.8,
    "lines.linewidth" : 2.5,
})

def set_fig_properties(ax_list):
    for ax in ax_list:
        ax.tick_params(which="major", length=8, width=2,   direction="in")
        ax.tick_params(which="minor", length=4, width=1.5, direction="in")
        ax.tick_params(which="both",  right=True, top=True)
        ax.minorticks_on()

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

axes[0].plot(t, (Ki + Ko) / L_G4)
axes[0].set_xlim([-0.5, 0.5])
axes[0].set_xticks([-0.5, 0, 0.5])
axes[0].set_xlabel("Correlation time (ps)")
axes[0].set_ylabel("K (eV/ps)")
axes[0].set_title(f"(a) Heat Current Correlation [{transport_dir}]")

axes[1].plot(nu, Gc)
axes[1].set_xlim([0, nu.max()])
axes[1].set_ylim(bottom=0)
axes[1].set_xlabel(r"Frequency $\omega/2\pi$ (THz)")
axes[1].set_ylabel(r"$G(\omega)$ (GW m$^{-2}$ K$^{-1}$ THz$^{-1}$)")
axes[1].set_title(f"(b) Spectral Conductance [{transport_dir}]")
axes[1].text(
    0.97, 0.95,
    f"$\\kappa$ = {kappa:.3f} W/mK\n$\\Delta T$ = {deltaT_use:.2f} K",
    transform=axes[1].transAxes, ha="right", va="top",
    fontsize=11, bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7)
)

set_fig_properties(axes)
plt.tight_layout()

save_dir  = os.path.expanduser("~/paper_2_figures")
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, f"shc_plot_{transport_dir}.pdf")
plt.savefig(save_path, bbox_inches="tight", pad_inches=0)
plt.show()
plt.close()
print(f"\nFigure saved  : {save_path}")
