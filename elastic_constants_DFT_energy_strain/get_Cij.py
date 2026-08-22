#!/usr/bin/env python3
"""
get_Cij.py
----------
Compute the elastic constants (C11, C12, C13, C33, C44, C66) of Bi2Te3
(or any rhombohedral/trigonal-symmetry crystal, point group -3m, with
the same strain sets) from the DFT energy-strain method. C66 is not
independent for this symmetry -- it is fixed by C66 = (C11-C12)/2 and
is reported for convenience, not fit from a separate strain set.

For each strain set (C11, C12, C13m, C33, C44), this fits
    (E(delta) - E0) / V0 = a * delta^2
and extracts the fit coefficient `a`. Each `a` is a known linear
combination of the true Cij (see COMBOS below, matching Table 3 of the
generation script). The combinations are then solved as a small linear
system to recover the individual Cij values.

Usage:
    python3 get_Cij.py /path/to/Bi2Te3

Expects this directory structure under the given root:
    <root>/C11/scf/delta_*/espresso.pwo
    <root>/C12/scf/delta_*/espresso.pwo
    <root>/C13m/scf/delta_*/espresso.pwo
    <root>/C33/scf/delta_*/espresso.pwo
    <root>/C44/scf/delta_*/espresso.pwo
"""
import sys
import os
import re
import glob
import numpy as np
from scipy.optimize import curve_fit

EV_A3_TO_GPA = 160.21766208
RY_TO_EV = 13.605693123
BOHR3_TO_ANG3 = 0.529177**3

# Strain-set -> combination of Cij that its quadratic fit coefficient
# corresponds to (a = coefficient of delta^2 in (E-E0)/V0).
# These match the comments in extract_elastic.py.
COMBOS = {
    "C11":  lambda C11, C12, C13, C33, C44: C11 + C12,
    "C12":  lambda C11, C12, C13, C33, C44: C11 - C12,
    "C13m": lambda C11, C12, C13, C33, C44: (C11 - 2*C13 + C33) / 2,
    "C33":  lambda C11, C12, C13, C33, C44: C33 / 2,
    "C44":  lambda C11, C12, C13, C33, C44: 2 * C44,
}


def read_final_energy(pwo_path):
    energy_ry = None
    with open(pwo_path, "r") as f:
        for line in f:
            if re.match(r"\s*!\s+total energy", line):
                energy_ry = float(line.split()[-2])
    if energy_ry is None:
        raise RuntimeError(f"No total energy found in {pwo_path}")
    return energy_ry * RY_TO_EV


def read_cell_volume(pwo_path):
    vol = None
    with open(pwo_path, "r") as f:
        for line in f:
            if "unit-cell volume" in line:
                m = re.search(r"=\s*([\d.]+)", line)
                if m:
                    vol = float(m.group(1)) * BOHR3_TO_ANG3
    if vol is None:
        raise RuntimeError(f"No unit-cell volume found in {pwo_path}")
    return vol


def fit_set(set_dir):
    """Return the quadratic fit coefficient `a` (eV/A^3) for one strain set."""
    delta_folders = sorted(glob.glob(os.path.join(set_dir, "delta_*")))
    if not delta_folders:
        raise FileNotFoundError(f"No delta_* folders in {set_dir}")

    records = []
    for d in delta_folders:
        pwo = os.path.join(d, "espresso.pwo")
        if not os.path.isfile(pwo):
            continue
        delta = float(os.path.basename(d).split("_")[1])
        E = read_final_energy(pwo)
        V = read_cell_volume(pwo)
        records.append((delta, E, V))

    if len(records) < 3:
        raise RuntimeError(f"Not enough completed runs in {set_dir}")

    records.sort(key=lambda x: x[0])
    records = np.array(records)
    deltas, energies, volumes = records[:, 0], records[:, 1], records[:, 2]

    zero_idx = np.argmin(np.abs(deltas))
    E0, V0 = energies[zero_idx], volumes[zero_idx]
    dE = (energies - E0) / V0

    def quadratic(x, a):
        return a * x**2

    popt, _ = curve_fit(quadratic, deltas, dE)
    return popt[0] * EV_A3_TO_GPA  # GPa


def solve_Cij(fit_values):
    """
    fit_values: dict with keys 'C11','C12','C13m','C33','C44' -> GPa
    Solves the linear system for the true C11, C12, C13, C33, C44.
    """
    # Unknowns: [C11, C12, C13, C33, C44]
    # Equations (rows match COMBOS, in the same order):
    #   C11 + C12                 = fit_values['C11']
    #   C11 - C12                 = fit_values['C12']
    #   (C11 - 2*C13 + C33)/2     = fit_values['C13m']
    #   C33/2                     = fit_values['C33']
    #   2*C44                     = fit_values['C44']
    A = np.array([
        [1,  1,  0, 0, 0],
        [1, -1,  0, 0, 0],
        [0.5, 0, -1, 0.5, 0],
        [0,  0,  0, 0.5, 0],
        [0,  0,  0, 0, 2],
    ])
    b = np.array([
        fit_values["C11"],
        fit_values["C12"],
        fit_values["C13m"],
        fit_values["C33"],
        fit_values["C44"],
    ])
    x = np.linalg.solve(A, b)
    C11, C12, C13, C33, C44 = x
    # C66 is not independent for trigonal/rhombohedral symmetry (R-3m,
    # point group -3m) -- it is fixed by C11 and C12, so no separate
    # strain set is needed for it.
    C66 = (C11 - C12) / 2
    return dict(C11=C11, C12=C12, C13=C13, C33=C33, C44=C44, C66=C66)


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 get_Cij.py /path/to/Bi2Te3")
        sys.exit(1)

    root = sys.argv[1]
    fit_values = {}
    print("Fitting each strain set...\n")
    for set_name in COMBOS:
        set_dir = os.path.join(root, set_name, "scf")
        try:
            a_gpa = fit_set(set_dir)
            fit_values[set_name] = a_gpa
            print(f"  {set_name:6s}  fit coefficient = {a_gpa:8.2f} GPa")
        except (FileNotFoundError, RuntimeError) as e:
            print(f"  {set_name:6s}  SKIPPED ({e})")

    missing = [s for s in COMBOS if s not in fit_values]
    if missing:
        print(f"\nCannot solve: missing completed sets {missing}")
        sys.exit(1)

    Cij = solve_Cij(fit_values)

    print("\n================= ELASTIC CONSTANTS (GPa) =================")
    for k in ["C11", "C12", "C13", "C33", "C44", "C66"]:
        note = "  (= (C11-C12)/2, not independently fit)" if k == "C66" else ""
        print(f"  {k} = {Cij[k]:8.2f} GPa{note}")
    print("=============================================================")


if __name__ == "__main__":
    main()
