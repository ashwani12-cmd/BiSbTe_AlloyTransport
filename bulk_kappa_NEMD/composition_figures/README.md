# κ∞ against composition — the three figures

Figure 5(a–c) of the manuscript, drawn from the extrapolated conductivities in
`../scaling_fits/` with the error bars of `../scaling_fits/kappa_inf_uncertainty.csv`.

| figure | what it shows |
|---|---|
| `fig_kappa_in_plane.pdf/.png` | κ∥ = (κ_X + κ_Y)/2 against Sb content, with a quadratic fit and an alloy-scattering (Klemens-type) fit |
| `fig_kappa_cross_plane.pdf/.png` | the same for κ⊥ = κ_Z |
| `fig_kappa_bulk_vs_literature.pdf/.png` | the bulk average (2κ∥ + κ⊥)/3 against five independent literature data sets |

| data file | content |
|---|---|
| `kappa_table_with_SE.csv` | κ∥, κ⊥ and the bulk average per composition, each with its standard error — the points plotted in all three figures |
| `literature_kappa_digitised.csv` | the comparison data of the third figure, digitised from the published sources; `series`, `kind` (points or curve), `x_Sb` as a fraction, `kappa` in W m⁻¹K⁻¹ |

The κ values are the thermostat-route intercepts; the error bars are `se_total` — the larger
of the block-averaged run-to-run scatter and the regression's own intercept SE. See
`../scaling_fits/README.md`.

## What the three figures say

Both components fall steeply from either end member into a broad minimum near equiatomic
composition — the signature of mass-disorder scattering. The measured drop over the first
20 % of substitution is steeper than either fitted form: at 20 % Sb the data point lies
below both the quadratic and the alloy-scattering curve, in-plane and cross-plane alike.

The two ends are not symmetric. Sb₂Te₃ conducts about 45 % better in-plane than Bi₂Te₃
(2.087 against 1.442 W m⁻¹K⁻¹) while being the more anisotropic of the two, so the minimum
sits off-centre and the two branches have different curvature.

Against experiment, the bulk average of this work follows the measured composition trend and
runs below Goldsmid and Rosi and above Birkholz through the alloy range, and below the
Boltzmann-transport calculation of Katcho *et al.* throughout. At the end members it is
below both measurement sets for Bi₂Te₃ (1.19 against 1.54 and 1.29) and between them for
Sb₂Te₃ (1.56 against 1.64 and 1.06).

## Literature sources in the third figure

- U. Birkholz, *Z. Naturforsch. A* **13**, 780 (1958). [10.1515/zna-1958-0910](https://doi.org/10.1515/zna-1958-0910)
- F. D. Rosi, B. Abeles, R. V. Jensen, *J. Phys. Chem. Solids* **10**, 191 (1959). [10.1016/0022-3697(59)90074-4](https://doi.org/10.1016/0022-3697(59)90074-4)
- H. J. Goldsmid, *J. Appl. Phys.* **32**, 2198 (1961). [10.1063/1.1777042](https://doi.org/10.1063/1.1777042)
- B. Poudel *et al.*, *Science* **320**, 634 (2008). [10.1126/science.1156446](https://doi.org/10.1126/science.1156446)
- N. A. Katcho, N. Mingo, D. A. Broido, *Phys. Rev. B* **85**, 115208 (2012). [10.1103/PhysRevB.85.115208](https://doi.org/10.1103/PhysRevB.85.115208)

Those points are **digitised from the published figures**, so they carry reading error on
top of whatever the original measurements carry; use them for the trend, not for
digit-by-digit comparison.
