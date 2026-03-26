import numpy as np
import matplotlib.pyplot as plt

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
    "legend.fontsize": 10,

    "axes.linewidth": 2.0,
    "lines.linewidth": 2.5,
    "lines.markersize": 8,

    "xtick.major.width": 2.0,
    "ytick.major.width": 2.0,
    "xtick.minor.width": 1.5,
    "ytick.minor.width": 1.5,
})

# ============================================================
# AXIS STYLE
# ============================================================
def apply_style(ax):
    ax.tick_params(axis='both', which='major',
                   direction='inout', length=8, width=2.0)
    ax.tick_params(axis='both', which='minor',
                   direction='in', length=6, width=1.5)
    ax.minorticks_on()

# ============================================================
# LOAD DATA
# ============================================================
data = np.loadtxt("loss.out", skiprows=1)

step = data[:, 0]
gen  = step / 100.0

total_loss = data[:, 1]
l1_loss    = data[:, 2]
l2_loss    = data[:, 3]

rmse_e_tr = data[:, 4]
rmse_f_tr = data[:, 5]
rmse_v_tr = data[:, 6]

# ============================================================
# PLOT: LOSS EVOLUTION (LEGEND INSIDE)
# ============================================================
fig, ax = plt.subplots(figsize=(5, 5))

ax.plot(gen, total_loss, lw=2, label="Total loss")
ax.plot(gen, l1_loss,    lw=2, label="L1 regularization")
ax.plot(gen, l2_loss,    lw=2, label="L2 regularization")

ax.plot(gen, rmse_e_tr,  lw=2, label="Energy (train)")
ax.plot(gen, rmse_f_tr,  lw=2, label="Force (train)")
ax.plot(gen, rmse_v_tr,  lw=2, label="Virial (train)")

ax.set_xscale("log")
ax.set_yscale("log")

apply_style(ax)

ax.set_xlabel("Generation / 100")
ax.set_ylabel("Loss functions")

# ============================================================
# LEGEND INSIDE (JOURNAL SAFE)
# ============================================================
leg = ax.legend(
    loc="best",
    frameon=True,
    fancybox=False,
    framealpha=0.85
)

leg.get_frame().set_edgecolor("black")
leg.get_frame().set_linewidth(1.0)


fig.tight_layout()
fig.savefig("loss_evolution.pdf", bbox_inches="tight", pad_inches=0)
plt.show()
plt.close(fig)

print("✅ Loss evolution plot saved (legend inside, publication-ready).")

