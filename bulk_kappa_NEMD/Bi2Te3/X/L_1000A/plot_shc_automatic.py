"""
plot_shc_auto.py
────────────────
Fully automated SHC plotter for GPUMD NEMD runs.
Reads ALL geometry/transport parameters from nemd_setup.txt — no hardcoding needed.
Works for X, Y, or Z transport directions and any composition.

Usage:
    python plot_shc_auto.py                    # looks for nemd_setup.txt in cwd
    python plot_shc_auto.py /path/to/run/dir   # pass a run directory explicitly
"""

import re
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from gpyumd.load import load_compute


# ============================================================
# 0. Resolve run directory
# ============================================================
run_dir = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
run_dir = os.path.abspath(run_dir)
print(f"Run directory : {run_dir}")

setup_file = os.path.join(run_dir, "nemd_setup.txt")
shc_file   = os.path.join(run_dir, "shc.out")

for f in (setup_file, shc_file):
    if not os.path.isfile(f):
        raise FileNotFoundError(f"Required file not found: {f}")


# ============================================================
# 1. Parse nemd_setup.txt
# ============================================================
def parse_setup(path):
    """Extract all needed scalars from nemd_setup.txt."""
    with open(path) as fh:
        text = fh.read()

    def grab(pattern, cast=float, flags=re.IGNORECASE):
        m = re.search(pattern, text, flags)
        if m is None:
            raise ValueError(f"Pattern not found in nemd_setup.txt: {pattern!r}")
        return cast(m.group(1))

    info = {}

    # Transport direction  (X / Y / Z)
    info["transport_dir"] = grab(r"Transport dir\s*:\s*([XYZ])", cast=str).upper()

    # Cell dimensions
    info["Lx"] = grab(r"Lx\s*:\s*([\d.]+)\s*A")
    info["Ly"] = grab(r"Ly\s*:\s*([\d.]+)\s*A")
    info["Lz"] = grab(r"Lz\s*:\s*([\d.]+)\s*A")

    # L_eff
    info["L_eff_A"] = grab(r"L_eff\s*:\s*([\d.]+)\s*A")

    # G4 length  — parse from group-boundaries table:
    #   GID 4   mid  G4   z_min   z_max   L(A) ...
    m4 = re.search(
        r"^\s*4\s+mid\s+G4\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
        text, re.MULTILINE | re.IGNORECASE
    )
    if m4 is None:
        raise ValueError("Cannot find GID 4 row in GROUP BOUNDARIES table.")
    info["L_G4"] = float(m4.group(3))   # L(A) column

    # Material / composition label (optional, used in title)
    m_mat = re.search(r"Material\s*:\s*(\S+)", text, re.IGNORECASE)
    info["material"] = m_mat.group(1) if m_mat else "Unknown"

    return info


params = parse_setup(setup_file)

transport_dir = params["transport_dir"]
Lx            = params["Lx"]
Ly            = params["Ly"]
Lz            = params["Lz"]
L_G4          = params["L_G4"]
L_eff_A       = params["L_eff_A"]
material      = params["material"]

# ── Cross-section area: the two dimensions perpendicular to transport ──────
cross_map = {
    "X": (Ly, Lz),
    "Y": (Lx, Lz),
    "Z": (Lx, Ly),
}
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
os.chdir(run_dir)   # load_compute / loadtxt look in cwd

data = np.loadtxt(shc_file)
corr = data[:499]
freq = data[499:]

t   = corr[:, 0]
Ki  = corr[:, 1]
Ko  = corr[:, 2]

nu  = freq[:, 0] / (2 * np.pi)
jwi = freq[:, 1]
jwo = freq[:, 2]


# ============================================================
# 3. Load temperature from compute.out
# ============================================================
compute  = load_compute(["temperature"])
T        = compute["temperature"]
ndata    = T.shape[0]
temp_ave = np.mean(T[int(ndata / 2) + 1:, 1:9], axis=0)
deltaT_actual = temp_ave[0] - temp_ave[-1]

print(f"temp_ave      : {np.round(temp_ave, 2)}")
print(f"deltaT        : {deltaT_actual:.4f} K")

deltaT_use = deltaT_actual if 5 < deltaT_actual < 100 else 20.0
if deltaT_use != deltaT_actual:
    print(f"  → deltaT out of (5,100) range; using fallback deltaT = {deltaT_use} K")


# ============================================================
# 4. Compute Gc and kappa
# ============================================================
Gc      = 1.6e4 * (jwi + jwo) / V / deltaT_use
G_total = np.trapz(Gc, nu)
L_eff_m = L_eff_A * 1e-10
kappa   = G_total * 1e9 * L_eff_m

print(f"Gc max        : {Gc.max():.6e} GW/m²/K/THz")
print(f"G_total       : {G_total * 1000:.4f} MW/m²/K")
print(f"kappa         : {kappa:.4f} W/mK")


# ============================================================
# 5. Plot
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

# Left panel — heat-current correlation
axes[0].plot(t, (Ki + Ko) / L_G4)
axes[0].set_xlim([-0.5, 0.5])
axes[0].set_xticks([-0.5, 0, 0.5])
axes[0].set_xlabel("Correlation time (ps)")
axes[0].set_ylabel("K (eV/ps)")
axes[0].set_title(f"(a) Heat Current Correlation [{transport_dir}]")

# Right panel — spectral conductance
axes[1].plot(nu, Gc)
axes[1].set_xlim([0, nu.max()])
axes[1].set_ylim(bottom=0)
axes[1].set_xlabel(r"Frequency $\omega/2\pi$ (THz)")
axes[1].set_ylabel(r"$G(\omega)$ (GW m$^{-2}$ K$^{-1}$ THz$^{-1}$)")
axes[1].set_title(f"(b) Spectral Conductance [{transport_dir}]")

set_fig_properties(axes)
plt.tight_layout()

# ── Save ──────────────────────────────────────────────────────────────────
save_dir  = os.path.expanduser("~/paper_2_figures")
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, f"shc_plot_{transport_dir}.pdf")
plt.savefig(save_path, bbox_inches="tight", pad_inches=0)
plt.show()
plt.close()
print(f"Figure saved  : {save_path}")
