#!/usr/bin/env python3
"""
make_fit_figures.py
--------------------
For each strain set (C11, C12, C13m, C33, C44), re-does the DFT
energy-strain quadratic fit and saves a labeled figure showing the
data points, the fit curve, and the fit equation + coefficient.

Usage:
    python3 make_fit_figures.py /path/to/Bi2Te3 [output_dir]
"""
import sys
import os
import re
import glob
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EV_A3_TO_GPA = 160.21766208
RY_TO_EV = 13.605693123
BOHR3_TO_ANG3 = 0.529177**3

SET_LABELS = {
    "C11":  r"$C_{11}+C_{12}$",
    "C12":  r"$C_{11}-C_{12}$",
    "C13m": r"$(C_{11}-2C_{13}+C_{33})/2$",
    "C33":  r"$C_{33}/2$",
    "C44":  r"$2\,C_{44}$",
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


def collect_set(set_dir):
    delta_folders = sorted(glob.glob(os.path.join(set_dir, "delta_*")))
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
    return np.array(records)


def fit_and_plot(set_name, set_dir, out_dir):
    records = collect_set(set_dir)
    deltas, energies, volumes = records[:, 0], records[:, 1], records[:, 2]

    zero_idx = np.argmin(np.abs(deltas))
    E0, V0 = energies[zero_idx], volumes[zero_idx]
    dE = (energies - E0) / V0  # eV/A^3

    def quadratic(x, a):
        return a * x**2

    popt, pcov = curve_fit(quadratic, deltas, dE)
    a = popt[0]
    a_err = np.sqrt(np.diag(pcov))[0]
    a_gpa = a * EV_A3_TO_GPA
    a_gpa_err = a_err * EV_A3_TO_GPA

    # R^2
    pred = quadratic(deltas, a)
    ss_res = np.sum((dE - pred) ** 2)
    ss_tot = np.sum((dE - dE.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

    # Plot
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 12.5,
        "axes.labelsize": 13.5,
        "axes.titlesize": 13,
        "legend.fontsize": 10.5,
    })
    fig, ax = plt.subplots(figsize=(5.6, 4.6))

    xfit = np.linspace(deltas.min() * 1.15, deltas.max() * 1.15, 300)
    yfit = quadratic(xfit, a) * 1000  # meV/A^3

    ax.scatter(deltas, dE * 1000, s=70, zorder=5, color="#c0392b",
               edgecolor="black", linewidth=0.8, label="DFT data")
    ax.plot(xfit, yfit, "--", color="#1a3a8f", lw=2, zorder=4, label="Quadratic fit")

    combo_label = SET_LABELS.get(set_name, set_name)
    eqn_text = (
        r"$\dfrac{E(\delta)-E_0}{V_0} = a\,\delta^2$" + "\n"
        rf"$a = {a*1000:.3f} \pm {a_err*1000:.3f}$ meV/$\mathrm{{\AA}}^3$" + "\n"
        rf"{combo_label} $= {a_gpa:.2f} \pm {a_gpa_err:.2f}$ GPa" + "\n"
        rf"$R^2 = {r2:.5f}$"
    )
    ax.text(0.03, 0.97, eqn_text, transform=ax.transAxes, fontsize=10,
            va="top", ha="left", fontfamily="serif", zorder=6,
            bbox=dict(boxstyle="round,pad=0.45", facecolor="white",
                      edgecolor="#888888", alpha=1.0))

    ax.set_xlabel(r"Strain $\delta$")
    ax.set_ylabel(r"$(E-E_0)/V_0$  (meV/$\mathrm{\AA}^3$)")
    ax.set_title(f"{set_name} strain set")
    ax.legend(loc="lower center", ncol=2, framealpha=0.95)
    ax.grid(alpha=0.3)
    ax.set_ylim(top=ax.get_ylim()[1] * 1.45)  # headroom so the box never covers data
    plt.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    out_pdf = os.path.join(out_dir, f"energy_strain_{set_name}_fit.pdf")
    out_png = os.path.join(out_dir, f"energy_strain_{set_name}_fit.png")
    plt.savefig(out_pdf, dpi=250, bbox_inches="tight")
    plt.savefig(out_png, dpi=250, bbox_inches="tight")
    plt.close()
    print(f"  {set_name}: a = {a_gpa:.2f} +/- {a_gpa_err:.2f} GPa, R2 = {r2:.5f}  -> {out_pdf}")
    return a_gpa


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 make_fit_figures.py /path/to/Bi2Te3 [output_dir]")
        sys.exit(1)

    root = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(root, "fit_figures")

    print("Generating fit figures for each strain set...\n")
    for set_name in ["C11", "C12", "C13m", "C33", "C44"]:
        set_dir = os.path.join(root, set_name, "scf")
        try:
            fit_and_plot(set_name, set_dir, out_dir)
        except (FileNotFoundError, RuntimeError) as e:
            print(f"  {set_name}: SKIPPED ({e})")

    print(f"\nAll figures saved to: {out_dir}")


if __name__ == "__main__":
    main()
