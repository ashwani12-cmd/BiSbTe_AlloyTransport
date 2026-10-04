#!/usr/bin/env python
"""
Figure: analytic interface conductance of Bi2Te3|Sb2Te3 vs temperature.

MODELS ONLY -- DMM (both Landauer conventions), AMM, phonon radiation limit.
No NEMD / GPUMD points: those belong in the separate NEMD figure.

Self-contained: the numbers below are hard-coded, so this runs anywhere with
just numpy + matplotlib.  They come from INTERFACE_CORRECT/s19_models.json
(script s19_mismatch_models.py), computed on the corrected, strain-balanced
R-3m interface cell at quasi-harmonically expanded lattice constants.

Style is the PUBLISHED paper style, lifted verbatim from the rcParams in your
plot_figures.ipynb -- serif/Times, inward ticks, 2.0 pt axes.  Don't swap it
for a seaborn/sans default, the figure has to sit next to Figs. 1-9.

    python plot_GK_models.py                 # writes fig_GK_models.{pdf,png}
    python plot_GK_models.py --out myname    # different basename
    python plot_GK_models.py --rad-both      # show the Bi2Te3-side limit too
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')          # drop this line if you want an interactive window
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- the data --
T = np.array([200, 300, 400, 500])                     # K

# all conductances in MW m^-2 K^-1
G_DMM_orig = np.array([33.96, 33.31, 32.47, 31.56])    # original Landauer convention
G_DMM_mod  = np.array([67.92, 66.61, 64.94, 63.13])    # modified = 2 x original
G_AMM      = np.array([92.45, 90.95, 89.03, 87.09])
G_rad_Sb   = np.array([95.75, 95.13, 93.50, 91.56])    # flux incident from Sb2Te3 side
G_rad_Bi   = np.array([105.85, 104.12, 101.87, 99.56]) # flux incident from Bi2Te3 side
# Table 'Radiation limit' column = G_rad_Sb, so that is the one plotted by default.

# flux-weighted transmission coefficients, for reference (not plotted)
alpha_DMM = np.array([0.321, 0.320, 0.319, 0.317])
alpha_AMM = np.array([0.873, 0.874, 0.874, 0.875])

# ------------------------------------------------- the published rc params --
PAPER_RC = {
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 16, "axes.labelsize": 18, "axes.titlesize": 18,
    "xtick.labelsize": 18, "ytick.labelsize": 18, "legend.fontsize": 14,
    "axes.linewidth": 2.0, "lines.linewidth": 2.5, "lines.markersize": 8,
    "xtick.major.width": 2.0, "ytick.major.width": 2.0,
    "xtick.minor.width": 1.5, "ytick.minor.width": 1.5,
}
_SIZE_KEYS = ("font.size", "axes.labelsize", "axes.titlesize",
              "xtick.labelsize", "ytick.labelsize", "legend.fontsize")
_WIDTH_KEYS = ("axes.linewidth", "lines.linewidth", "xtick.major.width",
               "ytick.major.width", "xtick.minor.width", "ytick.minor.width")


def use_paper_style(scale=1.0):
    """Reset to matplotlib defaults, then apply the paper style, scaled."""
    rc = dict(PAPER_RC)
    for k in _SIZE_KEYS:
        rc[k] = PAPER_RC[k] * scale
    for k in _WIDTH_KEYS:
        rc[k] = max(PAPER_RC[k] * scale, 0.5)
    rc["lines.markersize"] = max(PAPER_RC["lines.markersize"] * scale, 2.5)
    plt.rcParams.update(matplotlib.rcParamsDefault)
    plt.rcParams.update(rc)
    return scale


def apply_style(ax, scale=1.0, top=False, right=False, minor=True):
    """Inward major/minor ticks at the paper's weights."""
    ax.tick_params(axis='both', which='major', direction='in',
                   length=8 * scale, width=max(2.0 * scale, 0.5), top=top, right=right)
    ax.tick_params(axis='both', which='minor', direction='in',
                   length=6 * scale, width=max(1.5 * scale, 0.4), top=top, right=right)
    if minor:
        ax.minorticks_on()
    return ax


# ------------------------------------------------------------------- plot ---
ap = argparse.ArgumentParser()
ap.add_argument('--out', default='fig_GK_models', help='output basename')
ap.add_argument('--rad-both', action='store_true',
                help='also draw the Bi2Te3-side radiation limit')
ap.add_argument('--scale', type=float, default=1.0,
                help='scale all font/line sizes (e.g. 0.8 for a narrow column)')
a = ap.parse_args()

scale = use_paper_style(a.scale)
fig, ax = plt.subplots(figsize=(6.5, 5.5))

ax.plot(T, G_rad_Sb, ls=':', color='#7f7f7f', label='Radiation limit')
if a.rad_both:
    ax.plot(T, G_rad_Bi, ls=':', color='#c7c7c7',
            label=r'Radiation limit (Bi$_2$Te$_3$ side)')
ax.plot(T, G_AMM, ls='--', marker='s', color='#1f77b4', label='AMM')
ax.plot(T, G_DMM_mod, ls='-', marker='^', color='#2ca02c',
        label='DMM (mod. Landauer)')
ax.plot(T, G_DMM_orig, ls=(0, (1, 1)), marker='v', mfc='none', color='#2ca02c',
        label='DMM (original)')

ax.set_xlabel(r'Temperature (K)')
ax.set_ylabel(r'$G_{\mathrm{K}}$ (MW m$^{-2}$ K$^{-1}$)')
ax.set_xlim(170, 530)
ax.set_xticks(T)

lo = G_DMM_orig.min()
hi = (G_rad_Bi if a.rad_both else G_rad_Sb).max()
# headroom so the single-column legend sits above the curves, not on them
ax.set_ylim(lo - 8, hi + 55)
ax.legend(frameon=False, loc='upper left', ncol=1,
          handlelength=1.8, borderaxespad=0.4, labelspacing=0.3)
apply_style(ax, scale)

for ext in ('pdf', 'png'):
    fig.savefig(f'{a.out}.{ext}', dpi=300, bbox_inches='tight')
    print(f'wrote {a.out}.{ext}')
