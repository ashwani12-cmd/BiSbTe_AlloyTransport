"""
Plot R_K, dT_interface, and Finite-Size Correction (1/kappa vs 1/L)
as THREE SEPARATE figures (instead of one combined subplot).

Just edit the L, RK, dT, and inv_kappa_mKW arrays below with your
own values from analyze_interface.py output, then run:

    python3 plot_RK_three_figs.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# EDIT THESE WITH YOUR OWN VALUES (from analyze_interface.py output)
# ============================================================
L               = np.array([20.8, 31.0, 62.0])      # nm
RK              = np.array([7.40, 12.20, 16.30])    # x1e-9 m2K/W
dT              = np.array([3.36, 3.65, 3.17])       # K
invL            = 1.0 / L                            # nm^-1 (auto)
inv_kappa_mKW   = np.array([7.30, 9.55, 10.30])       # mK/W  (= 1/kappa, from your run)
# ============================================================

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 15.5,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "legend.fontsize": 11,
})

# ---------------------------------------------------------------
# Figure 1: Kapitza Resistance vs L
# ---------------------------------------------------------------
fig1, ax1 = plt.subplots(figsize=(5.5, 4.5))

ax1.plot(
    L, RK, 's-',
    color='purple',
    lw=3,
    ms=10,
    mfc='white',
    mew=2
)

mean_rk = RK.mean()

ax1.axhline(
    mean_rk,
    color='purple',
    ls='--',
    lw=1.5,
    label=f'mean={mean_rk:.3f} x1e-9 m2K/W'
)

# 7.40 — slightly to the RIGHT of the first marker
ax1.annotate(
    "7.40",
    (L[0], RK[0]),
    xytext=(15, 5),
    textcoords="offset points",
    ha="left",
    va="bottom",
    fontsize=10,
    fontweight="bold"
)

# 12.20 — above the middle marker
ax1.annotate(
    "12.20",
    (L[1], RK[1]),
    xytext=(0, 10),
    textcoords="offset points",
    ha="center",
    va="bottom",
    fontsize=10,
    fontweight="bold"
)

# 16.30 — INSIDE the plot, bottom-right
ax1.annotate(
    "16.30",
    (L[2], RK[2]),
    xytext=(-18, -20),
    textcoords="offset points",
    ha="left",
    va="top",
    fontsize=10,
    fontweight="bold"
)

ax1.set_xlabel('L [nm]',fontsize=18)
ax1.set_ylabel('R_K [x1e-9 m2K/W]',fontsize=18)

ax1.set_ylim(6.8, 16.8)

ax1.legend(fontsize=14)

plt.tight_layout()

plt.savefig("Fig1_RK_vs_L.png", dpi=300, bbox_inches="tight")
plt.savefig("Fig1_RK_vs_L.pdf", dpi=300, bbox_inches="tight")
plt.show()
plt.close()

print("Saved: Fig1_RK_vs_L.png / .pdf")


# ---------------------------------------------------------------
# Figure 2: Interface Temperature Jump vs L
# ---------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(5.5, 4.5))
ax2.plot(L, dT, '^-', color='darkorange', lw=2, ms=10, mfc='white', mew=2)
for l, d in zip(L, dT):
    ax2.annotate(f"{d:.2f}", (l, d), textcoords="offset points", xytext=(0, 10),
                 ha="center", fontsize=10, fontweight="bold")
ax2.set_xlabel('L (nm)')
ax2.set_ylabel('dT_interface (K)')
ax2.set_title('Interface Temperature Jump vs L')
ax2.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("Fig2_dT_vs_L.png", dpi=250)
plt.savefig("Fig2_dT_vs_L.pdf", dpi=250)
plt.close()
print("Saved: Fig2_dT_vs_L.png / .pdf")

# ---------------------------------------------------------------
# Figure 3: Finite-size correction (1/kappa vs 1/L)
# ---------------------------------------------------------------
cf = np.polyfit(invL, inv_kappa_mKW, 1)          # cf[0]=slope, cf[1]=intercept (mK/W)
pred = np.polyval(cf, invL)
ss_res = np.sum((inv_kappa_mKW - pred) ** 2)
ss_tot = np.sum((inv_kappa_mKW - inv_kappa_mKW.mean()) ** 2)
r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
intercept_mKW = cf[1]
kappa_inf_WmK = 1000.0 / intercept_mKW           # mK/W -> W/mK

fig3, ax3 = plt.subplots(figsize=(5.5, 4.5))
x_ext = np.linspace(0, invL.max() * 1.2, 300)
ax3.plot(invL, inv_kappa_mKW, 'o', color='C0', ms=11, mfc='white', mew=2,
         zorder=5, label='NEMD')
ax3.plot(x_ext, np.polyval(cf, x_ext), 'k--', lw=2, label=f'Fit R2={r2:.4f}')
ax3.plot(0, intercept_mKW, 'r*', ms=16, zorder=6,
         label=f'kappa_inf={kappa_inf_WmK:.3f} W/mK')
ax3.set_xlim(left=-0.002)
ax3.set_xlabel('1/L (nm-1)')
ax3.set_ylabel('1/kappa (mK/W)')
ax3.set_title('Finite-size correction')
ax3.legend()
plt.tight_layout()
plt.savefig("Fig3_FSC.png", dpi=250)
plt.savefig("Fig3_FSC.pdf", dpi=250)
plt.close()
print("Saved: Fig3_FSC.png / .pdf")

# ---------------------------------------------------------------
# Print summary table
# ---------------------------------------------------------------
print(f"\n{'L(nm)':>8} {'R_K(1e-9 m2K/W)':>18} {'dT_interface(K)':>17}")
for l, rk, d in zip(L, RK, dT):
    print(f"{l:>8.1f} {rk:>18.2f} {d:>17.2f}")
print(f"\nMean R_K = {mean_rk:.3f} x1e-9 m2K/W")
print(f"Fitted kappa_inf (1/L -> 0) = {kappa_inf_WmK:.3f} W/mK  (R2={r2:.4f})")
