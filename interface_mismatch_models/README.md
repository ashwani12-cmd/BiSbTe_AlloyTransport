# Interface mismatch models — AMM, DMM and the phonon radiation limit

Analytic interface thermal conductance of the Bi₂Te₃|Sb₂Te₃ junction, computed from the
NEP potential on a relaxed, strain-balanced cell at quasi-harmonically expanded geometries
for 200, 300, 400 and 500 K.

These are **lattice-dynamics calculations only** — no molecular dynamics enters this
folder. The NEMD interface runs live in `../interface_kappa_NEMD/`.

All energies and forces come from `../nep_train/nep.txt` (generation-100,000, md5
`545adcd43d8d8e19928a9f75937dee52`), symlinked here as `nep.txt`. If your checkout does not
preserve symlinks, copy the file in by hand.

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

## Results (`s19_models.json`), MW m⁻² K⁻¹

| T (K) | DMM (original) | DMM (mod. Landauer) | AMM | radiation limit (Sb side) | α_DMM | α_AMM |
|---|---|---|---|---|---|---|
| 200 | 33.96 | 67.92 | 92.45 | 95.75 | 0.321 | 0.873 |
| 300 | 33.31 | 66.61 | 90.95 | 95.13 | 0.320 | 0.874 |
| 400 | 32.47 | 64.94 | 89.03 | 93.50 | 0.319 | 0.874 |
| 500 | 31.56 | 63.13 | 87.09 | 91.56 | 0.317 | 0.875 |

Two conventions for the reference temperature drop are reported. The **original** form
refers G_K to the full reservoir difference Δ; the **modified Landauer** form uses the
local modal temperatures adjacent to the interface, for which ΔT_λ/Δ = 1/2 independently of
the transmission, so it is exactly 2× the original. `s19_mismatch_models.py` derives that
factor rather than hard-coding it, with an assertion, so the convention breaks loudly if
it is ever revised.

**α is referred to the Bi₂Te₃-side incident flux** (G_rad = 105.85, 104.12, 101.87,
99.56 MW m⁻² K⁻¹ at 200–500 K), whereas the radiation-limit column above is the *smaller*
Sb₂Te₃-side value, which is the binding bound. The two are different sides on purpose, so
α is not the ratio of the columns shown.

## Pipeline

Run in order, in this directory:

```bash
python s1_bulk.py             # free R-3m equilibrium of each compound
python s2_strain_balance.py   # the common in-plane lattice constant
python s5_phonon_check.py     # dynamical-stability check (no imaginary modes)
python s8_qha.py              # quasi-harmonic F(a,c,T) grid     [slowest stage]
python s19_mismatch_models.py # AMM + DMM + radiation limit  -> s19_models.json
python plot_GK_models.py      # -> fig_GK_models.pdf/.png
```

Requires `ase`, `calorine`, `phonopy`, `numpy`, `matplotlib`.

| file | what |
|---|---|
| `s1_bulk.py` → `s1_bulk_free.json`, `s1_*_free.xyz` | 15-atom R-3m cell, 3 quintuple layers, registry advancing on every plane **including the van der Waals gaps**. BFGS to fmax = 10⁻⁴ eV/Å. Free a = 4.389 (Bi₂Te₃) / 4.243 Å (Sb₂Te₃), misfit 3.4 % |
| `s2_strain_balance.py` → `s2_epitaxial.json`, `s2_*_epitaxial.xyz` | a scanned 4.22–4.42 Å in 0.01 Å steps, c and internals relaxed at fixed a, b. Balanced **a = 4.3162 Å**: −1.65 % on Bi₂Te₃, +1.73 % on Sb₂Te₃, residual σ_zz < 4×10⁻⁴ GPa |
| `s5_phonon_check.py` → `s5_phonons.json` | zero imaginary modes on either side; spectra end at 4.08 / 4.67 THz |
| `s8_qha.py` → `s8_qha_grid.json` | F(a,c,T) = E_static + F_vib on 5 in-plane × 6 out-of-plane points per compound, internals relaxed at each; F_vib from a 16×16×6 mesh. At each T the common a minimises the mean free energy and each side takes its own c(T) |
| `s19_mismatch_models.py` → `s19_models.json` | the three models. Force constants from 0.01 Å displacements in a 4×4×1 (240-atom) supercell, symmetrised; frequencies, eigenvectors and group velocities on a 24×24×8 mesh with symmetry reduction **off**; flux binned into 500 bins over 0–5.2 THz with modes below 0.05 THz dropped |
| `plot_GK_models.py` | the figure. Data hard-coded, so it runs standalone with only numpy + matplotlib |

## Two traps worth recording

Both were hit on the way to this version, and both give plausible-looking wrong answers:

1. **`primitive_matrix='auto'` puts phonopy in the 5-atom rhombohedral cell**, so `[0,0,q]`
   is *not* along Cartesian z and its period is not the hexagonal c. The q-vector must be
   built from the Cartesian k as `q_j = (a_j · k)/2π` over the primitive vectors, then
   `v = 2π·df/dk`. Getting this wrong gives sound speeds of 16–26 km/s instead of ~2 km/s.

2. **Sorting the three acoustic branches by frequency does not identify LA.** In Sb₂Te₃ the
   degenerate transverse pair lies *above* the longitudinal singleton, so "largest = LA"
   labels it backwards and the AMM then pairs LA on one side with TA on the other. The
   branches are classified by eigenvector instead — LA if the displacement is mostly along
   z. The mislabelling gives α_AMM = 0.931 instead of 0.873.

Sound velocities along z at 200 K (m/s): Bi₂Te₃ TA 2070.3 (×2), LA 2162.3; Sb₂Te₃
TA 2136.1 (×2), **LA 1884.5**. Sb₂Te₃ really does have LA slower than TA along z in this
potential (C₄₄ > C₃₃ for a strongly layered van der Waals solid), independently consistent
with the LAMMPS elastic constants in `../elastic_constants_LAMMPS/`. Acoustic-branch
linearity R² ≥ 0.999999.

## Scope and limitations

- **AMM's α carries no frequency dependence.** It is an elastic-continuum, long-wavelength
  quantity applied across the whole spectrum. That is the standard limitation of the model,
  not an extra approximation made here; the per-polarisation Γ_p are reported in
  `s19_models.json` under `AMM_per_pol` so the spread behind the average is visible.
- **The DMM is a lower bound for this pair, not an estimate.** Detailed balance caps α at
  0.5 even for two identical solids, and that cap bites hardest for closely matched
  materials — which Bi₂Te₃ and Sb₂Te₃ are, their acoustic impedances differing by only
  15–29 % (Z_Sb/Z_Bi = 0.845 transverse, 0.714 longitudinal).
- **No anharmonicity.** None of these models contains inelastic transmission, in which
  three- and four-phonon processes let energy cross the interface at a frequency different
  from the incident one. They are harmonic by construction.
- `s19_mismatch_models.py` also prints two reference columns (`NEMD_seed1`, `paper`) so the
  models can be compared against simulation at a glance. Those are **inputs for comparison
  only** — nothing in this folder computes them.
