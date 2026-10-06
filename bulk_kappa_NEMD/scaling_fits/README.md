# Finite-size scaling — the straight-line fits behind every κ∞

Every conductivity quoted for this system is an **extrapolation**, not a measurement: a
finite NEMD cell truncates the phonon mean free path, so κ(L) rises with system length and
the bulk value is recovered from

```
1/κ(L) = 1/κ∞ + Λ / (κ∞ L)
```

a straight line in (1/L, 1/κ) whose intercept is 1/κ∞ and whose slope-to-intercept ratio is
the effective mean free path Λ. This folder shows all 18 of those lines — 6 compositions ×
3 directions, 4 lengths each — so the extrapolation can be judged rather than taken on
trust.

| file | what |
|---|---|
| `Fig_scaling_fit_<comp>_<dir>.pdf/.png` | the 18 fits one at a time: the four NEMD points, the fit, and the intercept with its error bars |
| `Fig_scaling_fits.pdf/.png` | all 18 as one 6 × 3 panel figure, weak fits flagged in red |
| `scaling_fits.csv` | κ∞, SE(κ∞), relative SE, R², and Λ with its SE — for all three flux/geometry routes |
| `kappa_inf_uncertainty.csv` | the two uncertainties separately: block-averaged run-to-run (`se_stat`) and the regression's own (`se_regr`), with the larger one as `se_total` |
| `kappa_per_run.csv` | the 72 individual runs the fits are built from |
| `scaling_fits.py` | regenerates the CSVs and the panel figure from the raw data in `../<comp>/<dir>/L_*/` |

`python scaling_fits.py` (numpy + matplotlib) reproduces `scaling_fits.csv` exactly from the
gzipped run output in this repository. The per-series figures were drawn from the same fits
by the manuscript's revision-figure notebook; the numbers printed on them are the
`scaling_fits.csv` values.

## Two different error bars, and why both are drawn

Each fit has four points and therefore two degrees of freedom, so the regression's own
standard error on the intercept is a weak statement. It is drawn as the thin capped bar.
Behind it, the pale band is `se_total` — the larger of that and the block-averaged
run-to-run scatter from `kappa_inf_uncertainty.csv`. For the alloy cross-plane series the
run-to-run noise is 3–7× larger than the regression SE, so quoting the regression alone
would understate the uncertainty by that factor. **The tables quote `se_total`.**

## Four of the eighteen fits are weak, and they are not weak for the same reason

Flagged in red on the panel figure (R² < 0.90). The distinction matters more than the
number, because a low R² means different things on a flat line and on a steep one:

1. **20 % Sb, Z — R² = 0.779, but SE only 2.8 %.** The benign case. κ moves from 0.223 to
   0.247 W m⁻¹K⁻¹ across a fourfold length range, so the line is nearly flat (Λ = 3.4 nm)
   and there is little signal for R² to explain, while the intercept is tightly pinned. A
   low R² here says the series is already length-converged. **Quote the SE, not the R².**
2. **Bi₂Te₃, Z — R² = 0.824, SE 9.5 %.** All four runs are full length and pass the
   steady-state check (|imbalance| ≤ 0.084). The scatter is genuine run-to-run noise in a
   single realisation, not a defect: κ⊥(Bi₂Te₃) = 0.67 ± 0.06 W m⁻¹K⁻¹ is the honest result.
3. **80 % Sb, Y — R² = 0.888, SE 6.7 %.** All four runs clean. Mild curvature about the
   linear form; acceptable.
4. **80 % Sb, Z — R² = 0.899, SE 17.2 %.** The genuinely poor one. Two of its four runs are
   short (2850 and 3033 rows against a nominal 5000) and its 100 nm run has a source/sink
   power imbalance of −0.216, the largest in the whole 72-run set. Its fitted Λ = 33 nm is
   an outlier too (every other series gives 3–19 nm), the signature of an under-converged
   longest point dragging the slope up. κ⊥(80 % Sb) = 0.34 ± 0.06 W m⁻¹K⁻¹ carries real
   uncertainty and should be read with that caveat.

Three of the four are cross-plane, which is also the set the geometry convention moves most,
so the cross-plane row is the one that needs its error bars shown rather than bare numbers.

Every fit uses n = 4 lengths, so dof = 2 and the SEs are themselves uncertain. Eight of the
72 runs are shorter than the nominal 5000 rows, and the worst steady-state imbalances are
all cross-plane; both are visible per run in `kappa_per_run.csv`.

κ is the thermostat-power route with the true hexagonal cross-section throughout; see
`../README.md` for what that means and how it compares with the other two routes, which
`scaling_fits.csv` also carries.
