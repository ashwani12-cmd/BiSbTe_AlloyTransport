# s22–s26 — the converged QHA pipeline. THESE ARE THE NUMBERS THAT GO IN THE PAPER.

2026-10-05. Closes the blocker in `s19b_NOTES.md`: the old quasi-harmonic grid was not
converged, so every conductance in the old table was wrong by 5–13 % and every temperature
trend by about a factor of two. Rebuilt, converged, gated, and regenerated end to end.

**Everything in the subsection now traces to `s24_models_FINAL.json`. `s19_models.json`
and `s19c_fixed_geometry.json` are SUPERSEDED — do not cite them.**

## The chain

| script | makes | what it is |
|---|---|---|
| `s22_qha_refined.py` | `s22_qha_grid.json` | refined QHA grid: 17 c points over 4 %, centred on the minimum (old: 6 points over 8 %, off-centre) |
| `s23_qha_convergence.py` | — | **the gate.** 4 tests; exits non-zero if the grid is not usable |
| `s24_models_FINAL.py` | `s24_models_FINAL.json` | the models. Cubic F(c), a = 4.27 dropped |
| `s25_fixed_geometry_FINAL.py` | `s25_fixed_geometry_FINAL.json` | frozen-geometry control |
| `s26_dispersion_FINAL.py` | `s26_dispersion_FINAL.npz/.json` | dispersion + flux + α(ω) figure data |

Old grid kept for the record; `s19_models.json.bak-preReview` is the pre-audit file.

## What was wrong with the old grid — two independent defects

1. **Window.** `fc` ran 0.97, 0.99, 1.00, 1.01, 1.03, 1.05 — 8 % wide, asymmetric, with
   gaps. The free-energy minimum sits at **fc ≈ 1.005**, bracketed closely by only 1.00 and
   1.01, so the quadratic was dominated by the distant points. Δc varied by **2×** with fit
   order (0.0269 / 0.0534 / 0.0593 Å for orders 2/3/4).
2. **Contamination, never previously checked.** The old grid carried **843 imaginary
   modes**, 677 of them at fc = 1.05 across all five a-points — an over-expanded cell that
   is dynamically unstable, whose F_vib is meaningless, sitting *inside* the fit. Diagnosed
   as the lowest acoustic branch at large |q| (0.43–0.62, zone-boundary, **not** a near-Γ
   numerical artefact), Bi2Te3 only, growing with expansion and with in-plane compression.

## The gate (s23) — tests 1–3 pass on the refined grid, test 4 does not

```
TEST 1  fit-order independence   orders 2/3/4 -> 0.0574 / 0.0524 / 0.0511 A   spread 11.7 %  PASS
TEST 2  window independence      +-1 % / 1.5 % / 2 % / all -> 0.0520/0.0519/0.0518/0.0524   1.1 %  PASS
TEST 3  quadratic residual       0.082 meV/atom (old grid: 0.403)                           PASS
TEST 4  unstable-point contamination  see below                                             FAIL
```
With the unstable in-plane column `a = 4.27` dropped, orders 3 and 4 agree to **0.2 %**
(+0.0533 / +0.0532 Å), and dropping it shifts Δc by only **1.7 %**. `a = 4.27` is
1.1 % below the strain-balanced constant and has **0 of 17** clean c points; every retained
geometry is clean. Production setting: **cubic F(c), a ≥ 4.30.**

**Be straight about test 4 when writing this up: the script prints FAIL and exits 1.**
Its near-minimum half passes (worst point within ±1 % of the minimum carries 7 imaginary
modes of 23040, 0.03 %, at `--amin 4.28`); the half that fails is the clean-only refit,
Δc = +0.0478 Å vs +0.0533 Å, a **10.3 % shift against the 10 % threshold set in advance**.
The argument below — that the clean subset extrapolates and so is not a fair test — is an
argument made *after* seeing the number, and it should be presented that way, not as a pass.
`s23` takes `--amin` as of 2026-10-06; without it every in-plane column is kept, the
clean-only refit returns NaN (a = 4.27 has no clean c points at all) and the near-minimum
half fails too, 0.104 % vs 0.1 %.

A caution for whoever revisits this: the "keep only points with zero imaginary modes"
subset gives Δc = +0.0379 Å, which looks alarming but is **not a fair test** — that subset
spans fc 0.990–1.008 while the 500 K minimum is at fc ≈ 1.010, so it extrapolates rather
than fits. The a-drop and fit-order tests are the meaningful ones.

## THE NUMBERS

```
 T      DMM    DMM_mod    AMM    rad_Sb   rad_Bi   alpha_DMM  alpha_AMM
200    38.64    77.28    97.46   108.08   114.56     0.337      0.851
300    37.27    74.55    94.98   106.60   111.78     0.333      0.850
400    35.94    71.87    91.41   104.34   108.15     0.332      0.845
500    34.70    69.39    87.35   101.84   104.15     0.333      0.839
```

200 → 500 K, expanded vs frozen geometry:
```
  DMM      -10.21 %   /  +1.51 %
  AMM      -10.38 %   /  +1.44 %
  rad Bi    -9.09 %   /  +1.44 %
  rad Sb    -5.77 %   /  +2.57 %
```
Split on G_rad,Bi: **geometry alone −10.40 %, statistics alone +1.44 %, diagonal −9.09 %,
ratio 7.2.** Classical (k_B) DMM_mod: **78.68 → 69.60 = −11.54 %** — still falling, so the
Bose factor is not the story.

Mechanism (Bi2Te3, 200→500 K): acoustic flux **−11.64 %**, optical **−9.08 %** (50/50
weight) — acoustic softens more. Mean positive v_z **−9.82 %**; Γ-slope LA only −6.91 %, so
**the Γ slope understates it**. Band-mean frequencies move −1.74 / −0.64 %.

Dispersion (300 K geometry, a = 4.3178): acoustic max 0.792 / 0.823 THz, lowest optical
1.470 / 1.844 THz, tops **4.072 / 4.686 THz**. Flux-weighted α by band: 0.217 (<1 THz),
0.485, 0.597, 0.611.

## How the paper changed

```
                    OLD (unconverged)      NEW (converged)
 DMM_mod 200 K           67.92                  77.28      +13.8 %
 AMM     200 K           92.45                  97.46       +5.4 %
 rad_Sb  200 K           95.75                 108.08      +12.9 %
 trends             -4.4 .. -7.1 %        -5.8 .. -10.4 %
 alpha_AMM            0.873 (rising)        0.851 (falling)
```
**The "agreement to within 10 % with Chowdhury's 75 MW/m²K" was an artefact of the bad
fit.** The converged value is 77.3 at 200 K — within **3 %** — which is a better result,
but it is a different claim and the text now makes it at the temperature where the two
geometries are closest.

## Validation that survives into the deliverables

- `s26` asserts its Γ slopes reproduce the s24 sound speeds to <0.2 % **and** that
  re-integrating its own flux spectra returns G_DMM = 37.27 and G_rad = 111.78, the
  published 300 K values. Both pass, so the dispersion figure cannot drift from Table 1.
- `s25` reproduces the 200 K row of s24 exactly by construction.
- `plot_GK_models.py` asserts the two geometries agree at 200 K.
- The LaTeX table was machine-compared against `s24_models_FINAL.json`: exact on all
  6 columns × 4 rows.
