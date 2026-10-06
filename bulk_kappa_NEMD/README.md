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

**The figures are in `composition_figures/`** — κ∥, κ⊥ and the bulk average against Sb
content, with error bars, fitted curves and a comparison against five literature data sets,
plus the CSVs they are drawn from.

**The extrapolations themselves are in `scaling_fits/`** — all 18 straight-line fits of
1/κ against 1/L, one figure per composition and direction plus a 6 × 3 overview, with κ∞,
its standard error, R² and the effective mean free path in `scaling_fits.csv`. Four of the
eighteen are weak and that folder's README says which, and why each is weak for a different
reason.

Both components show a clear **alloy minimum near equiatomic composition** — κ∥ bottoms out
at x_Sb = 40 %, a factor of 2.0 below Bi₂Te₃ and 2.9 below Sb₂Te₃, which is mass-disorder
scattering doing what it should. The two end members are *not* equivalent: Sb₂Te₃ conducts
~45 % better in-plane than Bi₂Te₃ but is more anisotropic.

## Flux and geometry conventions

Two conventions enter a NEMD conductivity and both are worth stating explicitly, because
they move the answer by tens of per cent:

1. **Bin width.** The flux per bin must be divided by that bin's own volume. There are six
   mid bins between the source and the sink; a width of `L*(1-2*0.02-2*0.18)/8` spreads the
   same mid-region over eight, giving a width 4/3 too large.
2. **Cross-section.** `A_cross` as printed in `nemd_setup.txt` is `|b × c|`, the area of the
   parallelogram spanned by the two non-transport cell vectors. For a hexagonal cell
   (γ = 120°) the cross-section perpendicular to the transport axis is `V/L`, smaller by
   2/√3 = 1.1547. The true value is recomputed here from the NPT-equilibrated cell.

`kappa_composition.py` reports three conductivities per run so both conventions are visible
side by side:

| column | flux | geometry |
|---|---|---|
| `k_pub` | local `jp` | as published — both errors left in |
| `k_jp` | local `jp` | each bin's own width, true cross-section from the NPT cell |
| `k_th` | energy actually added at the source and removed at the sink | true cross-section; **no bin volume enters at all** |

**`k_th` is the route to quote**, and it is the one tabulated above: no bin volume enters it
at all. `k_pub` is retained because it is how the three routes can be compared on identical
runs — on these runs it returns 1.411 / 0.891 for Bi₂Te₃ and 2.004 / 0.605 for Sb₂Te₃.

In-plane, `k_pub` and `k_th` agree to within a few per cent, the two factors largely
cancelling in X and Y; the alloy minimum near equiatomic composition is the same either way.
Cross-plane carries the full 4/3, so κ⊥ is ~25 % lower on the `k_th` route and the
anisotropy κ∥/κ⊥ comes out at 2.1–4.2. The `paper_par` / `paper_perp` columns in the
script's last table and in `kappa_vs_composition.csv` hold the values as first published, so
any of these comparisons can be made directly from the CSV.

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
