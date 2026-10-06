#!/usr/bin/env python
"""
Figure: cross-plane phonon dispersion of the two sides of the Bi2Te3|Sb2Te3 junction,
with the two quantities that actually enter the mismatch models.

This is the figure Chowdhury et al. (ACS AMI 13, 4636 (2021)) show as their Fig. 2a --
the bulk dispersions along the transport direction -- on the CORRECTED R-3m
strain-balanced cell, and extended with the flux spectra and the DMM transmission.

Three panels sharing one frequency axis:
  (a) Gamma--Z dispersion, MIRRORED: Bi2Te3 opening to the left, Sb2Te3 to the right,
      Gamma common at the centre.  Gamma--Z is the transport direction (cross-plane).
  (b) the one-sided phonon flux Phi(omega) of each side -- Eq. (flux) in the text.
  (c) the DMM transmission alpha(omega) = Phi_2/(Phi_1+Phi_2) -- Eq. (dmm).

Reading it left to right: the two dispersions are nearly identical below ~1.5 THz
(hence the high AMM transmission), but Sb2Te3's spectrum runs 0.6 THz higher, and
every Bi2Te3 mode above 4.07 THz has nothing to scatter into.  alpha falls off a
cliff there.  That asymmetry is the whole reason alpha_DMM ~ 0.32 rather than 0.5.

Data: s26_dispersion_FINAL.npz (in this folder), written by s26_dispersion_FINAL.py
on the 300 K quasi-harmonic geometry (a = 4.3178 A) of the REFINED s22 grid, i.e. the
300 K row of s24_models_FINAL.json.
That script asserts the Gamma slopes reproduce the s19 sound speeds and that these
spectra re-integrate to the published G_DMM and G_rad.  Both pass.

Style is the PUBLISHED paper style -- same block as plot_GK_models.py.

    python plot_dispersion.py                # writes fig_dispersion.{pdf,png}
    python plot_dispersion.py --npz PATH     # point at the .npz elsewhere
"""
import argparse, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument('--npz', default='s26_dispersion_FINAL.npz')
ap.add_argument('--out', default='fig_dispersion')
ap.add_argument('--scale', type=float, default=1.0)
a = ap.parse_args()

ap_smooth = 5          # running-mean width, in 10.4 GHz bins -> ~52 GHz
d = np.load(a.npz)
fgrid, alpha_raw = d['fgrid'], d['alpha']
P1_raw, P2_raw = d['Bi_phi'], d['Sb_phi']
C_BI, C_SB = '#1f77b4', '#d62728'

def smooth(y, w=ap_smooth):
    """Running mean, for DISPLAY ONLY.

    The flux spectra are histograms of a 24x24x8 mesh into 500 bins, so they carry
    visible binning noise.  Every NUMBER quoted in the text and in s19_models.json is
    computed from the UNSMOOTHED spectra -- this touches the curves, nothing else.
    Edge bins are averaged over the shorter window rather than zero-padded, which
    would pull the ends towards zero."""
    if w <= 1: return y.copy()
    out = np.empty_like(y, dtype=float)
    h = w // 2
    for i in range(len(y)):
        lo, hi = max(0, i-h), min(len(y), i+h+1)
        out[i] = y[lo:hi].mean()
    return out

P1, P2 = smooth(P1_raw), smooth(P2_raw)
# alpha is recomputed from the smoothed fluxes so that panels (b) and (c) agree
den = P1 + P2
alpha = np.where(den > 0, P2/np.where(den > 0, den, 1.0), np.nan)

# Above Bi2Te3's cutoff there is NO incident flux, so alpha = Phi2/Phi2 = 1 is a 0/0
# artefact, not a transmission.  Mask it and shade the region instead.
f_cut = float(fgrid[P1_raw > 1e-3*P1_raw.max()][-1])
alpha_plot = np.where(P1_raw > 1e-3*P1_raw.max(), alpha, np.nan)

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
_SIZE = ("font.size","axes.labelsize","axes.titlesize","xtick.labelsize",
         "ytick.labelsize","legend.fontsize")
_WID  = ("axes.linewidth","lines.linewidth","xtick.major.width",
         "ytick.major.width","xtick.minor.width","ytick.minor.width")
rc = dict(PAPER_RC)
for k in _SIZE: rc[k] = PAPER_RC[k]*a.scale
for k in _WID:  rc[k] = max(PAPER_RC[k]*a.scale, 0.5)
plt.rcParams.update(matplotlib.rcParamsDefault); plt.rcParams.update(rc)

def style(ax, top=False, right=False, minor=True):
    ax.tick_params(axis='both', which='major', direction='in',
                   length=8*a.scale, width=max(2.0*a.scale,0.5), top=top, right=right)
    ax.tick_params(axis='both', which='minor', direction='in',
                   length=6*a.scale, width=max(1.5*a.scale,0.4), top=top, right=right)
    if minor: ax.minorticks_on()

FMAX = 5.0
fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.6),
                         gridspec_kw={'width_ratios': [1.55, 1.0, 1.0]})
axd, axf, axa = axes

# ---------------- (a) mirrored Gamma--Z dispersion --------------------------
# x runs -1..+1: Bi2Te3's Gamma->Z mapped onto -1..0 (reversed, so Gamma is at the
# centre), Sb2Te3's onto 0..+1.  Both axes are the SAME physical path; mirroring is
# a presentation choice that puts the two sides of the junction either side of Gamma.
for tag, col, sgn, nm in (('Bi', C_BI, -1.0, r'Bi$_2$Te$_3$'),
                          ('Sb', C_SB, +1.0, r'Sb$_2$Te$_3$')):
    k, bands = d[f'{tag}_k'], d[f'{tag}_bands']
    x = sgn * k / k.max()
    for b in range(bands.shape[1]):
        axd.plot(x, bands[:, b], color=col, lw=max(1.6*a.scale, 0.4),
                 label=nm if b == 0 else None)
axd.axvline(0.0, color='#999999', lw=max(1.2*a.scale, 0.4), zorder=0)
axd.set_xlim(-1, 1); axd.set_ylim(0, FMAX)
axd.set_xticks([-1, 0, 1]); axd.set_xticklabels([r'$Z$', r'$\Gamma$', r'$Z$'])
axd.set_ylabel(r'Frequency (THz)')
axd.set_xlabel(r'Cross-plane wavevector')
axd.legend(frameon=False, loc='upper center', ncol=2, handlelength=1.4,
           columnspacing=1.2, borderaxespad=0.3)
style(axd, top=True, right=True, minor=False)
axd.set_yticks(np.arange(0, FMAX+0.1, 1.0)); axd.minorticks_on()

# ---------------- (b) one-sided flux spectra --------------------------------
sc = 1e-18                       # Phi peaks at ~7.5e18 1/(m^2 s Hz)
axf.axhspan(f_cut, FMAX, color='#f0f0f0', zorder=0)
axf.plot(P1*sc, fgrid, color=C_BI, label=r'Bi$_2$Te$_3$ ($\Phi_1$)')
axf.plot(P2*sc, fgrid, color=C_SB, label=r'Sb$_2$Te$_3$ ($\Phi_2$)')
axf.set_ylim(0, FMAX); axf.set_xlim(0, 8.4)
axf.set_xticks([0, 2, 4, 6, 8])
axf.set_xlabel(r'$\Phi$ ($10^{18}$ m$^{-2}$s$^{-1}$Hz$^{-1}$)')
axf.legend(frameon=False, loc='upper right', handlelength=1.4, borderaxespad=0.3)
axf.set_yticklabels([])
style(axf, top=True, right=True)

# ---------------- (c) DMM transmission --------------------------------------
axa.axhspan(f_cut, FMAX, color='#f0f0f0', zorder=0)
axa.plot(alpha_plot, fgrid, color='#2ca02c')
axa.axvline(0.5, ls=':', color='#999999', lw=max(1.5*a.scale, 0.4))
axa.text(0.52, 4.78, r'no incident', fontsize=PAPER_RC['legend.fontsize']*a.scale*0.8,
         color='#666666', ha='center', va='center')
axa.text(0.52, 4.52, r'flux', fontsize=PAPER_RC['legend.fontsize']*a.scale*0.8,
         color='#666666', ha='center', va='center')
axa.set_ylim(0, FMAX); axa.set_xlim(0, 1.0)
axa.set_xticks([0, 0.5, 1.0])
axa.set_xlabel(r'$\alpha^{\mathrm{DMM}}_{1\rightarrow2}(\omega)$')
axa.set_yticklabels([])
style(axa, top=True, right=True)

for ax, tag in zip(axes, ('(a)', '(b)', '(c)')):
    ax.text(0.0, 1.015, tag, transform=ax.transAxes, ha='left', va='bottom',
            fontsize=PAPER_RC['axes.labelsize']*a.scale)
fig.subplots_adjust(wspace=0.12)

for ext in ('pdf', 'png'):
    fig.savefig(f'{a.out}.{ext}', dpi=300, bbox_inches='tight')
    print(f'wrote {a.out}.{ext}')

bi, sb = d['Bi_bands'].max(), d['Sb_bands'].max()
print(f"spectrum tops on the Gamma-Z path: Bi2Te3 {bi:.3f} THz, Sb2Te3 {sb:.3f} THz")
print(f"Bi2Te3 flux cutoff (shaded above): {f_cut:.3f} THz")
for lo, hi in ((0.05, 1.0), (0.05, 2.0)):
    m = (fgrid > lo) & (fgrid < hi)
    aw = P2_raw[m].sum()/(P1_raw[m]+P2_raw[m]).sum()     # from the RAW spectra
    print(f"flux-weighted alpha, {lo}-{hi} THz = {aw:.3f}")
print("  (Chowdhury et al. report 0.40-0.45 below 1 THz; ours is lower because this NEP")
print("   gives Sb2Te3 a much slower cross-plane LA, 1906 vs 2238 m/s -- see notes)")
