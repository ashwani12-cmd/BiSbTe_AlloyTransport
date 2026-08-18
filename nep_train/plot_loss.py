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

# ============================================================
# LIMIT TO 100k STEPS
# ============================================================
max_step = 300000
mask = step <= max_step

step = step[mask]
gen  = step / 100.0

total_loss = data[:, 1][mask]
l1_loss    = data[:, 2][mask]
l2_loss    = data[:, 3][mask]

rmse_e_tr = data[:, 4][mask]
rmse_f_tr = data[:, 5][mask]
rmse_v_tr = data[:, 6][mask]

# ============================================================
# PLOT
# ============================================================
fig, ax = plt.subplots(figsize=(5, 5))

ax.plot(gen, total_loss, label="Total loss")
ax.plot(gen, l1_loss,    label="L1 regularization")
ax.plot(gen, l2_loss,    label="L2 regularization")

ax.plot(gen, rmse_e_tr,  label="Energy (train)")
ax.plot(gen, rmse_f_tr,  label="Force (train)")
ax.plot(gen, rmse_v_tr,  label="Virial (train)")

# log scales
ax.set_xscale("log")
ax.set_yscale("log")

# apply style
apply_style(ax)

# labels
ax.set_xlabel("Generation / 100")
ax.set_ylabel("Loss functions")

# limit x-axis tightly to 100k
ax.set_xlim(gen.min(), gen.max())

# ============================================================
# LEGEND
# ============================================================
leg = ax.legend(
    loc="best",
    frameon=True,
    fancybox=False,
    framealpha=0.85
)

leg.get_frame().set_edgecolor("black")
leg.get_frame().set_linewidth(1.0)

# ============================================================
# SAVE
# ============================================================
fig.tight_layout()
fig.savefig("loss_evolution_100k.pdf", bbox_inches="tight", pad_inches=0)

plt.show()
plt.close(fig)

print("✅ Loss evolution plot (0–100k steps) saved successfully.")
