# Bulk thermal conductivity by NEMD — the full composition series at 300 K

Non-equilibrium MD on (Bi₁₋ₓSbₓ)₂Te₃ for **x_Sb = 0, 20, 40, 60, 80 and 100 %**, each along
the three crystallographic directions **X, Y, Z**, each at **four system lengths**
(25, 50, 75, 100 nm) for finite-size extrapolation. 72 runs, 6 × 3 × 4, complete.

```
<composition>/<direction>/L_<length>A/
```

X and Y are in-plane (a- and b-axis), Z is cross-plane (c-axis, across the van der Waals
gaps). The end members Bi₂Te₃ (x = 0) and Sb₂Te₃ (x = 100) are included in the same grid, so
the series is continuous from one end member to the other.

Potential: `../nep_train/nep.txt`, generation-100,000, md5
`545adcd43d8d8e19928a9f75937dee52`. All runs are 1 ns NPT equilibration followed by 5 ns of
Langevin-thermostat NEMD at T = 300 K, ΔT = 20 K (source 320 K / sink 280 K), on a
cross-section targeted at 120 nm² for every direction and composition.

## Results at 300 K (W m⁻¹ K⁻¹)

Reproduce with `python kappa_composition.py` (numpy only):

| x_Sb (%) | κ∥ (in-plane) | κ⊥ (cross-plane) | anisotropy κ∥/κ⊥ |
|---|---|---|---|
| 0 (Bi₂Te₃) | 1.442 | 0.673 | 2.14 |
| 20 | 0.859 | 0.250 | 3.44 |
| 40 | 0.712 | 0.223 | 3.19 |
| 60 | 0.752 | 0.241 | 3.12 |
| 80 | 0.972 | 0.339 | 2.87 |
| 100 (Sb₂Te₃) | 2.087 | 0.503 | 4.15 |

κ∥ = (κ_X + κ_Y)/2, κ⊥ = κ_Z, each from the 1/κ vs 1/L intercept over the four lengths.
Per-run values, gradients, areas and fit quality are in `kappa_composition.csv`; the series
above is `kappa_vs_composition.csv`.

Both components show a clear **alloy minimum near equiatomic composition** — κ∥ bottoms out
at x_Sb = 40 %, a factor of 2.0 below Bi₂Te₃ and 2.9 below Sb₂Te₃, which is mass-disorder
scattering doing what it should. The two end members are *not* equivalent: Sb₂Te₃ conducts
~45 % better in-plane than Bi₂Te₃ but is more anisotropic.

## Read this before using the numbers: the published table used a different recipe

The conductivities printed in the paper were computed from the local virial heat flux `jp`
with two geometry errors, and the numbers above are **not** the published ones:

1. **Bin width.** The published post-processing divided the per-bin flux by a bin volume
   built as `L*(1-2*0.02-2*0.18)/8`, i.e. it split the whole mid-region into 8 bins when
   there are only 6 mid bins between source and sink. The bin width is 4/3 too large.
2. **Cross-section.** `A_cross` as printed in `nemd_setup.txt` is `|b × c|`, the area of the
   parallelogram spanned by the two non-transport cell vectors. For a hexagonal cell
   (γ = 120°) the true cross-section perpendicular to the transport axis is `V/L`, smaller by
   2/√3 = 1.1547.

`kappa_composition.py` reports three conductivities per run so the correction is auditable
rather than asserted:

| column | flux | geometry |
|---|---|---|
| `k_pub` | local `jp` | as published — both errors left in |
| `k_jp` | local `jp` | each bin's own width, true cross-section from the NPT cell |
| `k_th` | energy actually added at the source and removed at the sink | true cross-section; **no bin volume enters at all** |

**Quote `k_th`.** It is the route that does not depend on a bin volume, and it is the one
tabulated above. `k_pub` is kept because it is the only way to check the diagnosis: running
the published recipe on these same runs returns the published numbers (Bi₂Te₃ 1.411 / 0.891
against 1.37 / 0.88 in print, Sb₂Te₃ 2.004 / 0.605 against 2.02 / 0.69), which is what
identifies the two factors above as the whole of the discrepancy.

What this changes:

- **In-plane κ survives.** `k_pub` and `k_th` agree to a few per cent in X and Y, because the
  4/3 and the area factor partly cancel there. The alloy minimum and every conclusion drawn
  from the in-plane series stand.
- **Cross-plane κ must be restated.** Z carries the full 4/3, so κ⊥ drops by ~25 %. The
  anisotropy κ∥/κ⊥ is **2.1–4.2**, not the 1.5–2.0 quoted in the paper. In particular
  κ⊥(Bi₂Te₃) = 0.673, so the claimed agreement with Witting *et al.* (~0.85) does not hold.

The `paper_par` / `paper_perp` columns in the script's last table and in
`kappa_vs_composition.csv` carry the published values for exactly this comparison. Note also
that Table 5 of the paper labels the composition axis "% Bi" where it means **% Sb**.

## What is in each run directory

| file | what |
|---|---|
| `nemd_setup.txt` | the full geometry: supercell, L_transport, A_cross **as printed** (see above), the 8 group boundaries and bin centres. The analysis takes the bin grid from here |
| `compute.out.gz` | the NEMD data: per-group temperature, the three `jp` blocks, and the accumulated source/sink thermostat energies in the last two columns. 5000 rows = 5 ns at 1 ps |
| `thermo.out.gz` | thermodynamic output including the instantaneous cell vectors, from which the **true** cross-section is recomputed after NPT |
| `run.in`, `submit.sh` | the GPUMD input and the job script, verbatim |
| `stoichiometry.txt` | alloy runs only: target and achieved Sb fraction, atom counts, and the **substitution random seed** |
| `shc.out` | spectral heat current (GPUMD `compute_shc`), where the run produced one |
| `out.dat` | per-run post-processing summary as it was generated at the time |

`compute.out` and `thermo.out` are gzipped to keep the repository a sensible size; numpy
reads `.gz` directly, so nothing needs unpacking and the analysis script takes either form.
The Bi₂Te₃ directories predate this and still hold their files uncompressed, along with the
original per-direction `post_processing_nemd_Jp.py` and its figures.

**`model.xyz` and `relaxed.xyz` are deliberately not included** — 60 structures of 85,000+
atoms is 2 GB. They are regenerable instead: `build_nemd_cells_AS_USED.py` is the generator
exactly as it was run, it reads the DFT-relaxed conventional cells in
`../BiSbTe_alloy_endpoints/` (`espresso_{Bi2Te3,Sb2Te3}_conventional.pwi`, which reproduce
the supercell dimensions in `nemd_setup.txt` to within 2×10⁻⁴ relative), and the Bi→Sb
substitution seed is `RANDOM_SEED = 1`, recorded per run in `stoichiometry.txt`. Note that
this generator is also where the `A_cross` of point 2 above comes from — its `get_A_cross()`
returns `|b × c|`. It is shipped unmodified, because it is what produced these runs.
