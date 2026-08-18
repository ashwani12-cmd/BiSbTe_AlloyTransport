import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# ============================================================
# JOURNAL-QUALITY GLOBAL SETTINGS
# ============================================================
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",

    "font.size": 16,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 18,

    "axes.linewidth": 2.0,
    "lines.linewidth": 2.5,
    "lines.markersize": 8,

    "xtick.major.width": 2.0,
    "ytick.major.width": 2.0,
    "xtick.minor.width": 1.5,
    "ytick.minor.width": 1.5,
})

def apply_style(ax):
    ax.tick_params(axis='both', which='major',
                   direction='in', length=8, width=2.0,
                   top=False, right=False)
    ax.tick_params(axis='both', which='minor',
                   direction='in', length=6, width=1.5,
                   top=False, right=False)
    ax.minorticks_on()

def parity_axes(ax, x, y, nticks=5, pad=0.05):
    vmin = min(x.min(), y.min())
    vmax = max(x.max(), y.max())
    dv = vmax - vmin
    vmin -= pad * dv
    vmax += pad * dv
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)
    ax.plot([vmin, vmax], [vmin, vmax], '--', color='gray', lw=2)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=nticks))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=nticks))
    ax.set_aspect("auto")

# ============================================================
# LOAD TEST DATA
# ============================================================
energy = np.loadtxt("energy_test.out")
force  = np.loadtxt("force_test.out")
stress = np.loadtxt("stress_test.out")

# ENERGY
dft_e = energy[:, 0]
nep_e = energy[:, 1]
rmse_e = np.sqrt(np.mean((dft_e - nep_e) ** 2))

# FORCE
nep_f = force[:, 0:3]
dft_f = force[:, 3:6]
rmse_f = np.sqrt(np.mean((dft_f - nep_f) ** 2))
labels_f = ["Fx", "Fy", "Fz"]
colors_f = ["tab:blue", "tab:orange", "tab:green"]

# STRESS
nep_s = stress[:, 0:6]
dft_s = stress[:, 6:12]
mask = np.all(dft_s > -1e5, axis=1)
nep_s = nep_s[mask]; dft_s = dft_s[mask]
rmse_s = np.sqrt(np.mean((dft_s - nep_s) ** 2))
labels_s = ["xx", "yy", "zz", "xy", "yz", "zx"]
colors_s = ["tab:blue", "tab:orange", "tab:green",
            "tab:red", "tab:purple", "tab:brown"]

print(f"Test  Energy RMSE: {rmse_e*1000:.3f} meV/atom")
print(f"Test  Force  RMSE: {rmse_f*1000:.3f} meV/Å")
print(f"Test  Stress RMSE: {rmse_s:.4f} GPa")

# ============================================================
# 1 — ENERGY PARITY (TEST)
# ============================================================
fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(dft_e, nep_e, s=18, label="energy")
parity_axes(ax, dft_e, nep_e)
apply_style(ax)
ax.set_xlabel("DFT energy [eV/atom]")
ax.set_ylabel("NEP energy [eV/atom]")
ax.text(0.18, 0.10, f"RMSE = {rmse_e*1000:.2f} meV/atom",
        transform=ax.transAxes)
leg = ax.legend(loc="upper left", frameon=True, fancybox=False, framealpha=0.85)
leg.get_frame().set_edgecolor("black"); leg.get_frame().set_linewidth(1.0)
fig.tight_layout()
fig.savefig("parity_energy_test_100k.pdf", bbox_inches="tight", pad_inches=0)
fig.savefig("parity_energy_test_100k.png", dpi=300, bbox_inches="tight", pad_inches=0)
plt.close(fig)
print("Saved: parity_energy_test_100k")

# ============================================================
# 2 — FORCE PARITY (TEST)
# ============================================================
fig, ax = plt.subplots(figsize=(5, 5))
for i in range(3):
    ax.scatter(dft_f[:, i], nep_f[:, i], s=10, color=colors_f[i], label=labels_f[i])
parity_axes(ax, dft_f.flatten(), nep_f.flatten())
apply_style(ax)
ax.set_xlabel("DFT force [eV/Å]")
ax.set_ylabel("NEP force [eV/Å]")
ax.text(0.20, 0.10, f"RMSE = {rmse_f*1000:.2f} meV/Å",
        transform=ax.transAxes)
leg = ax.legend(loc="upper left", frameon=True, fancybox=False, framealpha=0.85)
leg.get_frame().set_edgecolor("black"); leg.get_frame().set_linewidth(1.0)
fig.tight_layout()
fig.savefig("parity_force_test_100k.pdf", bbox_inches="tight", pad_inches=0)
fig.savefig("parity_force_test_100k.png", dpi=300, bbox_inches="tight", pad_inches=0)
plt.close(fig)
print("Saved: parity_force_test_100k")

# ============================================================
# 3 — STRESS PARITY (TEST)
# ============================================================
fig, ax = plt.subplots(figsize=(5, 5))
for i in range(6):
    ax.scatter(dft_s[:, i], nep_s[:, i], s=12, color=colors_s[i], label=labels_s[i])
parity_axes(ax, dft_s.flatten(), nep_s.flatten())
apply_style(ax)
ax.set_xlabel("DFT stress [GPa]")
ax.set_ylabel("NEP stress [GPa]")
ax.text(0.20, 0.10, f"RMSE = {rmse_s:.3f} GPa",
        transform=ax.transAxes)
leg = ax.legend(loc="upper left", frameon=True, fancybox=False, framealpha=0.85)
leg.get_frame().set_edgecolor("black"); leg.get_frame().set_linewidth(1.0)
fig.tight_layout()
fig.savefig("parity_stress_test_100k.pdf", bbox_inches="tight", pad_inches=0)
fig.savefig("parity_stress_test_100k.png", dpi=300, bbox_inches="tight", pad_inches=0)
plt.close(fig)
print("Saved: parity_stress_test_100k")

print("\nAll test parity plots done!")
