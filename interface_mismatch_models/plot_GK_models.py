#!/usr/bin/env python
"""
Figure: analytic interface conductance of Bi2Te3|Sb2Te3 vs temperature.

MODELS ONLY -- DMM (both Landauer conventions), AMM, phonon radiation limit.
No NEMD / GPUMD points: those belong in the separate NEMD figure.

TWO PANELS (new 2026-10-05):
  (a) absolute G_K(T) on the quasi-harmonically expanded geometries.
  (b) the SAME models normalised to their 200 K value, computed two ways --
      with the quasi-harmonic geometry (what panel (a) and the table show) and
      with the lattice FROZEN at its 200 K configuration.  Panel (b) is the
      whole point: frozen, the models rise and saturate, which is the textbook
      harmonic result; expanded, they fall.  The fall is the geometry, and
      without this panel a referee has to take that on trust.

Self-contained: the numbers below are hard-coded, so this runs anywhere with
just numpy + matplotlib.  Sources:
  panel (a) + QHA curves  -> s24_models_FINAL.json  (s24_models_FINAL.py)
  frozen-geometry curves  -> s25_fixed_geometry_FINAL.json
both on the corrected, strain-balanced R-3m cell, on the REFINED QHA grid (s22) with the
dynamically unstable a = 4.27 column dropped and a cubic F(c) fit -- see s19b_NOTES.md.
The earlier numbers (DMM_mod 67.92 etc.) came from an unconverged quadratic fit on the
old s8 grid and are superseded.

Style is the PUBLISHED paper style, lifted verbatim from the rcParams in your
plot_figures.ipynb -- serif/Times, inward ticks, 2.0 pt axes.  Don't swap it
for a seaborn/sans default, the figure has to sit next to Figs. 1-9.

    python plot_GK_models.py                 # writes fig_GK_models.{pdf,png}
    python plot_GK_models.py --out myname    # different basename
    python plot_GK_models.py --rad-both      # show the Bi2Te3-side limit too in (a)
    python plot_GK_models.py --single        # panel (a) only, the old one-panel figure
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')          # drop this line if you want an interactive window
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- the data --
T = np.array([200, 300, 400, 500])                     # K

# --- quasi-harmonically expanded geometry (s19_models.json) -----------------
# all conductances in MW m^-2 K^-1
G_DMM_orig = np.array([38.64, 37.27, 35.94, 34.70])    # original Landauer convention
G_DMM_mod  = np.array([77.28, 74.55, 71.87, 69.39])    # modified = 2 x original
G_AMM      = np.array([97.46, 94.98, 91.41, 87.35])
G_rad_Sb   = np.array([108.08, 106.60, 104.34, 101.84])# flux incident from Sb2Te3 side
G_rad_Bi   = np.array([114.56, 111.78, 108.15, 104.15])# flux incident from Bi2Te3 side
# Table 'Radiation limit' column = G_rad_Sb, so that is the one plotted by default.

# --- lattice FROZEN at the 200 K geometry (s19c_fixed_geometry.json) --------
# a = 4.3172 A, cQL(Bi2Te3) = 10.3086 A, cQL(Sb2Te3) = 9.8520 A at every T;
# only dn/dT varies.  The 200 K entries are identical to the QHA ones above by
# construction -- that identity is the control's own validation.
F_DMM_orig = np.array([38.64, 39.03, 39.16, 39.23])
F_DMM_mod  = np.array([77.28, 78.05, 78.32, 78.45])
F_AMM      = np.array([97.46, 98.38, 98.71, 98.86])
F_rad_Bi   = np.array([114.56, 115.64, 116.03, 116.21])

# flux-weighted transmission coefficients, for reference (not plotted).
# alpha_AMM RISES with T under expansion: the AMM factor 4Z1Z2/(Z1+Z2)^2 depends
# on the impedance RATIO, which barely moves when both sides soften together.
alpha_DMM = np.array([0.3373, 0.3334, 0.3324, 0.3331])
alpha_AMM = np.array([0.8507, 0.8496, 0.8452, 0.8387])
alpha_AMM_frozen = 0.8507                              # exactly constant, same phonons

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

C_AMM, C_DMM, C_RAD = '#1f77b4', '#2ca02c', '#7f7f7f'


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


def panel_absolute(ax, scale, rad_both=False):
    """(a) G_K(T), quasi-harmonically expanded geometry."""
    ax.plot(T, G_rad_Sb, ls=':', color=C_RAD, label='Radiation limit')
    if rad_both:
        ax.plot(T, G_rad_Bi, ls=':', color='#c7c7c7',
                label=r'Radiation limit (Bi$_2$Te$_3$ side)')
    ax.plot(T, G_AMM, ls='--', marker='s', color=C_AMM, label='AMM')
    ax.plot(T, G_DMM_mod, ls='-', marker='^', color=C_DMM,
            label='DMM (mod. Landauer)')
    ax.plot(T, G_DMM_orig, ls=(0, (1, 1)), marker='v', mfc='none', color=C_DMM,
            label='DMM (original)')
    ax.set_xlabel(r'Temperature (K)')
    ax.set_ylabel(r'$G_{\mathrm{K}}$ (MW m$^{-2}$ K$^{-1}$)')
    ax.set_xlim(170, 530)
    ax.set_xticks(T)
    lo = G_DMM_orig.min()
    hi = (G_rad_Bi if rad_both else G_rad_Sb).max()
    # headroom so the single-column legend sits above the curves, not on them
    ax.set_ylim(lo - 8, hi + 46)
    ax.legend(frameon=False, loc='upper left', ncol=1,
              handlelength=1.8, borderaxespad=0.4, labelspacing=0.3)
    apply_style(ax, scale)


def panel_normalised(ax, scale):
    """(b) the same models / their 200 K value: expanded vs frozen geometry.

    DMM's two conventions differ by an exact factor 2, so they normalise onto
    one curve and only one DMM line is drawn.
    """
    n = lambda y: y / y[0]
    ax.axhline(1.0, color='#cccccc', lw=max(1.2 * scale, 0.4), zorder=0)
    # frozen geometry -- open markers, dashed
    ax.plot(T, n(F_AMM), ls='--', marker='s', mfc='none', color=C_AMM,
            label='AMM, fixed geometry')
    ax.plot(T, n(F_DMM_mod), ls='--', marker='^', mfc='none', color=C_DMM,
            label='DMM, fixed geometry')
    # quasi-harmonic -- filled markers, solid
    ax.plot(T, n(G_AMM), ls='-', marker='s', color=C_AMM,
            label='AMM, expanded')
    ax.plot(T, n(G_DMM_mod), ls='-', marker='^', color=C_DMM,
            label='DMM, expanded')
    ax.set_xlabel(r'Temperature (K)')
    ax.set_ylabel(r'$G_{\mathrm{K}}(T)\,/\,G_{\mathrm{K}}(200\,\mathrm{K})$')
    ax.set_xlim(170, 530)
    ax.set_xticks(T)
    # derived from the data, never hard-coded: the trends changed by ~2x once the QHA
    # grid was converged, and a fixed limit silently clipped the 500 K points
    ys = np.concatenate([n(F_AMM), n(F_DMM_mod), n(G_AMM), n(G_DMM_mod)])
    pad = 0.12 * (ys.max() - ys.min())
    ax.set_ylim(ys.min() - pad, ys.max() + 2.4 * pad)
    ax.legend(frameon=False, loc='lower left', ncol=1,
              handlelength=1.8, borderaxespad=0.4, labelspacing=0.3)
    apply_style(ax, scale)


# ------------------------------------------------------------------- plot ---
ap = argparse.ArgumentParser()
ap.add_argument('--out', default='fig_GK_models', help='output basename')
ap.add_argument('--rad-both', action='store_true',
                help='also draw the Bi2Te3-side radiation limit in panel (a)')
ap.add_argument('--single', action='store_true',
                help='panel (a) only (the previous one-panel figure)')
ap.add_argument('--scale', type=float, default=1.0,
                help='scale all font/line sizes (e.g. 0.8 for a narrow column)')
a = ap.parse_args()

# the control's validating identity, asserted so a bad edit cannot pass silently
for nm, qha, frz in (('DMM', G_DMM_orig, F_DMM_orig), ('AMM', G_AMM, F_AMM),
                     ('rad', G_rad_Bi, F_rad_Bi)):
    assert abs(qha[0] - frz[0]) < 0.02, \
        f"{nm}: the 200 K point must be identical in both geometries, got {qha[0]} vs {frz[0]}"

scale = use_paper_style(a.scale)

if a.single:
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    panel_absolute(ax, scale, a.rad_both)
else:
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.5))
    panel_absolute(axes[0], scale, a.rad_both)
    panel_normalised(axes[1], scale)
    # tags ABOVE the axes: inside the frame they collide with panel (a)'s legend
    for ax, tag in zip(axes, ('(a)', '(b)')):
        ax.text(0.0, 1.015, tag, transform=ax.transAxes, ha='left', va='bottom',
                fontsize=PAPER_RC['axes.labelsize'] * a.scale)
    fig.subplots_adjust(wspace=0.30)

for ext in ('pdf', 'png'):
    fig.savefig(f'{a.out}.{ext}', dpi=300, bbox_inches='tight')
    print(f'wrote {a.out}.{ext}')

print("\n200 -> 500 K relative change (%):")
for nm, qha, frz in (('DMM', G_DMM_orig, F_DMM_orig), ('AMM', G_AMM, F_AMM),
                     ('radiation, Bi side', G_rad_Bi, F_rad_Bi)):
    print(f"  {nm:20s} expanded {100*(qha[-1]/qha[0]-1):+6.2f}    "
          f"fixed geometry {100*(frz[-1]/frz[0]-1):+6.2f}")
