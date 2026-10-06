# Interface mismatch models — AMM, DMM and the phonon radiation limit

Analytic interface thermal conductance of the Bi₂Te₃|Sb₂Te₃ junction, computed from the
NEP potential on a relaxed, strain-balanced R-3m cell at quasi-harmonically expanded
geometries for 200, 300, 400 and 500 K.

These are **lattice-dynamics calculations only** — no molecular dynamics enters this
folder. The NEMD interface runs live in `../interface_kappa_NEMD/`.

All energies and forces come from `../nep_train/nep.txt` (generation-100,000, md5
`545adcd43d8d8e19928a9f75937dee52`), symlinked here as `nep.txt`. If your checkout does not
preserve symlinks, copy the file in by hand.

> **This folder was regenerated on 2026-10-05.** The quasi-harmonic grid behind the first
> version was not converged: every conductance was 5–14 % low and every temperature trend
> about a factor of two too weak. The numbers below come from the rebuilt grid (`s22`) and
> the re-run models (`s24`). The superseded `s19_*` files have been removed rather than left
> to be cited by mistake; what was wrong with them is recorded in `s19b_NOTES.md` and
> `s24_FINAL_NOTES.md`, and the old grid itself (`s8_qha_grid.json`) is kept because the
> convergence gate compares against it.

## What is computed

All three models share the Landauer form

```
G_K = ∫ ħω Φ₁(ω) α₁→₂(ω) (∂n/∂T) dω
```

and differ only in the transmission coefficient α:

| model | α | note |
|---|---|---|
| phonon radiation limit | 1 | rigorous upper bound on any harmonic model |
| AMM (acoustic mismatch) | 4Z₁Z₂cosθ₁cosθ₂ / (Z₁cosθ₁+Z₂cosθ₂)², angle-averaged | specular; impedance Z = ρv; frequency-independent |
| DMM (diffuse mismatch) | Φ₂/(Φ₁+Φ₂) | detailed balance; frequency-resolved, no angles |

Φ(ω) is the one-sided phonon flux taken from the **full NEP dispersion**, not a Debye
approximation.

## Results (`s24_models_FINAL.json`), MW m⁻² K⁻¹

| T (K) | DMM (original) | DMM (mod. Landauer) | AMM | radiation limit (Sb side) | α_DMM | α_AMM |
|---|---|---|---|---|---|---|
| 200 | 38.64 | 77.28 | 97.46 | 108.08 | 0.337 | 0.851 |
| 300 | 37.27 | 74.55 | 94.98 | 106.60 | 0.333 | 0.850 |
| 400 | 35.94 | 71.87 | 91.41 | 104.34 | 0.332 | 0.845 |
| 500 | 34.70 | 69.39 | 87.35 | 101.84 | 0.333 | 0.839 |

Two conventions for the reference temperature drop are reported. The **original** form
refers G_K to the full reservoir difference Δ; the **modified Landauer** form uses the local
modal temperatures adjacent to the interface, for which ΔT_λ/Δ = 1/2 independently of the
transmission, so it is exactly 2× the original. `s24_models_FINAL.py` derives that factor
rather than hard-coding it, with an assertion, so the convention breaks loudly if it is ever
revised.

**α is referred to the Bi₂Te₃-side incident flux** (G_rad = 114.56, 111.78, 108.15,
104.15 MW m⁻² K⁻¹ at 200–500 K), whereas the radiation-limit column above is the *smaller*
Sb₂Te₃-side value, which is the binding bound. The two are different sides on purpose, so α
is not the ratio of the columns shown. Classical (k_B) counterparts of every column are in
the JSON under `*_classical`.

## All three models fall with temperature, and that is geometric

Over 200 → 500 K: DMM −10.2 %, AMM −10.4 %, radiation limit −9.1 % (Bi side) / −5.8 %
(Sb side). Harmonic models at *fixed* geometry are expected to rise and then saturate, so
the direction needs evidence rather than assertion. `s25_fixed_geometry_FINAL.py` is that
control: it freezes the lattice at its 200 K quasi-harmonic solution and varies only ∂n/∂T.

| 200 → 500 K | expanded (`s24`) | frozen (`s25`) |
|---|---|---|
| DMM | −10.21 % | +1.51 % |
| AMM | −10.38 % | +1.44 % |
| radiation limit, Bi side | −9.09 % | +1.44 % |

Frozen, the models rise and saturate, which is the textbook result (and what Reddy *et al.*,
Appl. Phys. Lett. **87**, 211908 (2005) show for the DMM on exact Born–von Kármán
dispersions). Expanded, they fall. The geometry term beats the occupation term by a factor
of about 7 (geometry alone −10.40 %, statistics alone +1.44 %, on G_rad,Bi).

**It is the flux, not the transmission.** α_AMM moves 0.851 → 0.839 and α_DMM is flat to
1 %; the loss is in Φ(ω), because every group velocity drops as the van der Waals gaps
widen. Decomposed for Bi₂Te₃: acoustic flux −11.6 %, optical −9.1 % at roughly equal weight,
mean positive v_z −9.8 %. **The Γ-slope sound speeds fall only 6.9 % and therefore understate
the effect** — quoting them alone cannot account for the total. Panel (b) of
`fig_GK_models.pdf` is this comparison; on an absolute axis a 10 % effect is invisible.

## Pipeline

Run in order, in this directory:

```bash
python s1_bulk.py                    # free R-3m equilibrium of each compound
python s2_strain_balance.py          # the common in-plane lattice constant
python s5_phonon_check.py            # dynamical-stability check (no imaginary modes)
python s22_qha_refined.py            # quasi-harmonic F(a,c,T) grid     [slowest stage]
python s23_qha_convergence.py --amin 4.28   # the gate: is that grid converged?
python s24_models_FINAL.py           # AMM + DMM + radiation limit -> s24_models_FINAL.json
python s25_fixed_geometry_FINAL.py   # frozen-geometry control
python s26_dispersion_FINAL.py       # dispersion + flux + α(ω) at the 300 K geometry
python plot_GK_models.py             # -> fig_GK_models.pdf/.png
python plot_dispersion.py            # -> fig_dispersion.pdf/.png
```

Requires `ase`, `calorine`, `phonopy`, `numpy`, `matplotlib`. `s8_qha.py` is the old grid,
kept only so the gate has something to compare against.

| file | what |
|---|---|
| `s1_bulk.py` → `s1_bulk_free.json`, `s1_*_free.xyz` | 15-atom R-3m cell, 3 quintuple layers, registry advancing on every plane **including the van der Waals gaps**. BFGS to fmax = 10⁻⁴ eV/Å. Free a = 4.389 (Bi₂Te₃) / 4.243 Å (Sb₂Te₃), misfit 3.4 % |
| `s2_strain_balance.py` → `s2_epitaxial.json`, `s2_*_epitaxial.xyz` | a scanned 4.22–4.42 Å in 0.01 Å steps, c and internals relaxed at fixed a, b. Balanced **a = 4.3162 Å**: −1.65 % on Bi₂Te₃, +1.73 % on Sb₂Te₃, residual σ_zz < 4×10⁻⁴ GPa |
| `s5_phonon_check.py` → `s5_phonons.json` | zero imaginary modes on either side |
| `s22_qha_refined.py` → `s22_qha_grid.json` | F(a,c,T) = E_static + F_vib, **17 c points over 4 %, centred on the minimum**, 5 in-plane points per compound, internals relaxed at each; F_vib from a 16×16×6 mesh |
| `s23_qha_convergence.py` | the convergence gate, 4 tests, exits non-zero if the grid is not usable. `--amin 4.28` reproduces the production setting; see the next section for what it does and does not establish |
| `s24_models_FINAL.py` → `s24_models_FINAL.json` | the three models. Cubic F(c), in-plane column a = 4.27 dropped. Force constants from 0.01 Å displacements in a 4×4×1 (240-atom) supercell, symmetrised; frequencies, eigenvectors and group velocities on a 24×24×8 mesh with symmetry reduction **off**; flux binned into 500 bins over 0–5.2 THz with modes below 0.05 THz dropped |
| `s25_fixed_geometry_FINAL.py` → `s25_fixed_geometry_FINAL.json` | the frozen-geometry control. Reproduces the 200 K row of `s24` exactly by construction, which is its own validation |
| `s26_dispersion_FINAL.py` → `s26_dispersion_FINAL.npz/.json` | Γ–Z dispersion, flux spectra and α(ω) at the 300 K geometry. Asserts its Γ slopes reproduce the `s24` sound speeds to < 0.2 % **and** that re-integrating its own spectra returns G_DMM = 37.27 and G_rad = 111.78, the tabulated 300 K values |
| `plot_GK_models.py`, `plot_dispersion.py` | the two figures. `plot_GK_models.py` has its data hard-coded, so it runs standalone with only numpy + matplotlib; it asserts the two geometries agree at 200 K |
| `AMM_DMM_EXPLAINED.md` | what the two models assume, where each is the right one, and the per-polarisation ingredients |
| `s19b_NOTES.md`, `s24_FINAL_NOTES.md` | the audit that found the unconverged grid, and the regeneration |

A few scripts still print the *old* output filenames at the end (`s24_models_FINAL.py` says
"wrote s19_models.json"). The scripts are shipped byte-identical to the ones that produced
the published numbers, so those messages have been left alone; the filename at the top of
each script is the one actually written.

## What the gate establishes, and what it does not

Run at the production setting, `s23_qha_convergence.py --amin 4.28`:

```
TEST 1  fit-order independence   orders 2/3/4 -> 0.0579 / 0.0533 / 0.0532 A  spread 8.5 %  PASS
TEST 2  window independence      +-1 / 1.5 / 2 % / all -> 0.0531/0.0534/0.0532/0.0533  0.5 %  PASS
TEST 3  quadratic residual       0.0799 meV/atom  (old grid: 0.4028)                    PASS
TEST 4  unstable-point contamination                                                    FAIL
```

Tests 1–3 are the point of the rebuild and they pass: Δc(200→500 K) no longer depends on
the fit (orders 3 and 4 agree to 0.2 %, against a factor of two on the old grid) or on where
the fitting window is cut.

**Test 4 fails, and the script says so rather than being tuned until it passes.** It has two
halves. The near-minimum half passes — the worst grid point within ±1 % of the minimum
carries 7 imaginary modes of 23040 (0.03 %). The half that fails refits using only grid
points with *zero* imaginary modes anywhere, which gives Δc = +0.0478 Å against +0.0533 Å,
a 10.3 % shift against a 10 % threshold. That subset spans fc = 0.990–1.008 while the 500 K
minimum sits at fc ≈ 1.010, so it extrapolates rather than fits, and we do not think it is a
fair test; the honest statement is that it is a criterion we set in advance, missed by
0.3 percentage points, and are arguing against after the fact. Readers should know that.

The test we do rely on for the unstable points is dropping the whole offending in-plane
column: a = 4.27 is 1.1 % below the strain-balanced constant and has **0 of 17** clean c
points, and removing it moves Δc by only **1.7 %** (0.0524 → 0.0533 Å). Run the gate with
and without `--amin 4.28` to reproduce that. Without the flag every column is kept, which is
the harder configuration, and test 4's first half then fails as well (0.104 % vs the 0.1 %
threshold).

## Two traps worth recording

Both were hit on the way here, and both give plausible-looking wrong answers:

1. **`primitive_matrix='auto'` puts phonopy in the 5-atom rhombohedral cell**, so `[0,0,q]`
   is *not* along Cartesian z and its period is not the hexagonal c. The q-vector must be
   built from the Cartesian k as `q_j = (a_j · k)/2π` over the primitive vectors, then
   `v = 2π·df/dk`. Getting this wrong gives sound speeds of 16–26 km/s instead of ~2 km/s.

2. **Sorting the three acoustic branches by frequency does not identify LA.** In Sb₂Te₃ the
   degenerate transverse pair lies *above* the longitudinal singleton, so "largest = LA"
   labels it backwards and the AMM then pairs LA on one side with TA on the other. The
   branches are classified by eigenvector instead — LA if the displacement is mostly along z.

Sound velocities along z at 200 K (m/s): Bi₂Te₃ TA 2145.0 (×2), LA 2287.9; Sb₂Te₃
TA 2242.6 (×2), **LA 1918.8**. Sb₂Te₃ really does have LA slower than TA along z in this
potential (C₄₄ > C₃₃ for a strongly layered van der Waals solid), independently consistent
with the LAMMPS elastic constants in `../elastic_constants_LAMMPS/`. Acoustic-branch
linearity R² ≥ 0.999999. Per-polarisation AMM at 200 K: LA Z = 18.283 / 12.562 MRayl,
Γ = 0.9709, **no critical angle** (Sb₂Te₃'s LA is the slower one); TA (×2) Z = 17.141 /
14.681 MRayl, Γ = 0.7906, θ_c = 73.0°, i.e. TA is the lossy channel. ρ = 7991.2 (Bi) /
6546.6 (Sb) kg m⁻³.

## Scope and limitations

- **AMM's α carries no frequency dependence.** It is an elastic-continuum, long-wavelength
  quantity applied across the whole spectrum. That is the standard limitation of the model,
  not an extra approximation made here; the per-polarisation Γ_p are in
  `s24_models_FINAL.json` under `AMM_per_pol` so the spread behind the average is visible.
- **The DMM is a lower bound for this pair, not an estimate.** Detailed balance caps α at
  0.5 even for two identical solids, and that cap bites hardest for closely matched
  materials — which Bi₂Te₃ and Sb₂Te₃ are, their acoustic impedances differing by only
  14–31 % (Z_Sb/Z_Bi = 0.857 transverse, 0.687 longitudinal).
- **No anharmonicity.** None of these models contains inelastic transmission, in which three-
  and four-phonon processes let energy cross the interface at a frequency different from the
  incident one. They are harmonic by construction.
- **The potential's weak spot is exactly the direction of transport here.** Against
  literature values the NEP gets Bi₂Te₃ C₃₃ almost exactly right but Sb₂Te₃ C₃₃ about 30 %
  too soft, so the Sb₂Te₃ cross-plane velocities — and through them α and every conductance
  above — inherit that error.
- `s24_models_FINAL.py` also prints two reference columns (`NEMD_seed1`, `paper`) so the
  models can be compared against simulation at a glance. Those are **inputs for comparison
  only** — nothing in this folder computes them.
