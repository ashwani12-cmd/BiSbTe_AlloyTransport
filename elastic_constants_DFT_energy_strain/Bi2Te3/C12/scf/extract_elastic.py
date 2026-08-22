#!/usr/bin/env python3
"""
extract_elastic.py
------------------
Run from:  C11/scf/   (or any set folder)
Reads espresso.pwo from each delta_*/  subfolder.
Extracts total energy and cell volume from QE output.
Fits E(delta) = E0 + a*delta^2 to get elastic constant.

Change SET_NAME and ELASTIC_LABEL per run.
"""

import numpy as np
import glob
import os
import re
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# ============================================================
# USER SETTINGS
# ============================================================
SET_NAME      = "C11"   # folder name — also used for plot title
ELASTIC_LABEL = "C11"   # label for the extracted constant

# Conversion
EV_A3_TO_GPA = 160.21766208

# ============================================================
# HELPERS: parse QE .pwo
# ============================================================
def read_final_energy(pwo_path):
    """Return final total energy in eV from espresso.pwo."""
    energy_ry = None
    with open(pwo_path, "r") as f:
        for line in f:
            # last occurrence of '!    total energy'
            if re.match(r"\s*!\s+total energy", line):
                energy_ry = float(line.split()[-2])   # in Ry
    if energy_ry is None:
        raise RuntimeError(f"No total energy found in {pwo_path}")
    return energy_ry * 13.605693123  # Ry -> eV

def read_cell_volume(pwo_path):
    """Return cell volume in Å³ from the last 'unit-cell volume' line."""
    vol = None
    with open(pwo_path, "r") as f:
        for line in f:
            if "unit-cell volume" in line:
                # format:  unit-cell volume          =     XXX.XXXX (a.u.)^3
                m = re.search(r"=\s*([\d.]+)", line)
                if m:
                    vol_au3 = float(m.group(1))   # in bohr^3
                    vol = vol_au3 * 0.529177**3   # bohr^3 -> Å³
    if vol is None:
        raise RuntimeError(f"No unit-cell volume found in {pwo_path}")
    return vol

# ============================================================
# STEP 1: collect delta, energy, volume
# ============================================================
delta_folders = sorted(glob.glob("delta_*"))
if not delta_folders:
    raise SystemExit("No delta_* folders found. Run from C11/scf/ (or the right set/scf dir).")

records = []
for d in delta_folders:
    pwo = os.path.join(d, "espresso.pwo")
    if not os.path.isfile(pwo):
        print(f"  SKIP {d} — espresso.pwo missing")
        continue
    delta = float(d.split("_")[1])   # e.g. delta_+0.002 -> +0.002
    try:
        E   = read_final_energy(pwo)
        vol = read_cell_volume(pwo)
        records.append((delta, E, vol))
        print(f"  delta={delta:+.4f}  E={E:.8f} eV  V={vol:.4f} Å³")
    except RuntimeError as e:
        print(f"  ERROR {d}: {e}")

records.sort(key=lambda x: x[0])
records = np.array(records)          # shape (N, 3)
deltas  = records[:, 0]
energies= records[:, 1]
volumes = records[:, 2]

# ============================================================
# STEP 2: E0 and V0 from delta = 0 folder
# ============================================================
zero_idx = np.argmin(np.abs(deltas))
E0 = energies[zero_idx]
V0 = volumes[zero_idx]
print(f"\nReference  delta={deltas[zero_idx]:+.4f}")
print(f"E0 = {E0:.8f} eV")
print(f"V0 = {V0:.4f} Å³")

# ============================================================
# STEP 3: energy density
# ============================================================
dE = (energies - E0) / V0   # eV/Å³

# ============================================================
# STEP 4: quadratic fit  dE = a * delta^2
# ============================================================
def quadratic(x, a):
    return a * x**2

popt, _ = curve_fit(quadratic, deltas, dE)
a = popt[0]

# ============================================================
# STEP 5: elastic constant  (see Table 3 for prefactor)
# ============================================================
# C11 set: d²E/dδ² / V0 = 2*(C11+C12)  -> prefactor in fit = (C11+C12)
# But for the pure C11 extraction, the fit gives:
#   dE = (C_combo) * delta^2   where C_combo depends on the strain set
# The prefactor mapping (from Table 3):
#   C11 set  -> a = C11+C12         => printed as "combo"
#   C12 set  -> a = C11-C12
#   C33 set  -> a = C33/2
#   C44 set  -> a = 2*C44
#   C13p set -> a = (C11+2C13+C33)/2
#   C13m set -> a = (C11-2C13+C33)/2
#   Bulk set -> a = (2C11+2C12+C33+4C13)/2
C_GPa = a * EV_A3_TO_GPA

# ============================================================
# OUTPUT
# ============================================================
print("\n================= RESULTS =================")
print(f"Fit coefficient  a  = {a:.6e} eV/Å³")
print(f"{ELASTIC_LABEL} combo = {C_GPa:.2f} GPa")
print("  (see Table 3 for exact C_ij combination)")
print("===========================================\n")

# ============================================================
# STEP 6: plot
# ============================================================
xfit = np.linspace(deltas.min(), deltas.max(), 300)
yfit = quadratic(xfit, a)

plt.figure(figsize=(5, 4))
plt.scatter(deltas, dE * 1000, s=60, zorder=5, label="DFT data")
plt.plot(xfit, yfit * 1000, "--", label=f"Fit: {ELASTIC_LABEL} combo = {C_GPa:.1f} GPa")
plt.xlabel("Strain δ")
plt.ylabel(r"$(E - E_0)\,/\,V_0$  (meV/Å³)")
plt.title(f"{SET_NAME} strain set")
plt.legend()
plt.grid(True)
plt.tight_layout()

out_pdf = f"energy_strain_{SET_NAME}.pdf"
plt.savefig(out_pdf, bbox_inches="tight")
print(f"Saved plot -> {out_pdf}")
plt.show()
