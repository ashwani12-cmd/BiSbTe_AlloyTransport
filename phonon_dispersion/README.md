# Phonon dispersion — DFT vs NEP

Data behind **Figure 3** of the manuscript: the phonon dispersion of Bi₂Te₃ and Sb₂Te₃
along Γ–M–K–Γ, computed with Quantum ESPRESSO and with the NEP potential in identical
cells.

## Layout

```
<compound>/
  band_dft.yaml      phonopy band output, DFT force constants
  band_nep.yaml      phonopy band output, NEP force constants
  band.conf          the phonopy band settings used for BOTH
  phonopy_disp.yaml  the displacement set
  POSCAR             unit cell (Bi2Te3)
  espresso*.pwi      the QE input
  band_comparison.png
  plot_comparison.py
  BiSbTe_alloy_dispersion_pipeline.sh
```

## Settings (`band.conf`, identical for DFT and NEP)

```
DIM  = 3 3 2
BAND = 0 0 0   0.5 0 0   0.3333 0.3333 0   0 0 0        # Γ → M → K → Γ
BAND_POINTS = 51
FC_SYMMETRY = .TRUE.
```

Both force sets are evaluated on **the same displaced supercells**, so the DFT and NEP
curves differ only through the forces.

### What is deliberately not included

The raw force files (`FORCE_SETS`, `FORCE_CONSTANTS`, per-displacement outputs) are **not**
published here. During preparation their provenance could not be established to our own
satisfaction: the accompanying `phonopy.yaml` does not parse, and the force magnitudes in
`FORCE_SETS` are roughly two orders of magnitude smaller than a 0.02 Å displacement should
produce, which points to a unit-convention or staleness problem in those particular files
rather than in the dispersions themselves. Rather than publish data we cannot vouch for, we
ship the band results, which **have** been verified (below), together with every input
needed to regenerate them. Anyone wishing to reproduce the force constants can do so from
`POSCAR`, `phonopy_disp.yaml` and `band.conf` with the potential in `../nep_train/`.

## ⚠ Two things that will confuse you if nobody says them

**1. `band_nep.yaml` is byte-identical to `qe_forces/band.yaml` in the source tree.**
That looks as though the "NEP" file is really the QE one. It is not. The directory is named
after the **Quantum ESPRESSO displacement structures**, not the force source; the
`band.yaml` sitting inside it was produced from NEP forces on those structures. This was
verified independently by recomputing the Bi₂Te₃ NEP dispersion from scratch with
`phonopy` + `calorine` on the same q-path: it reproduces `band_nep.yaml` to an RMS of
**0.012 THz** (maximum 3.806 vs 3.807 THz) and differs from `band_dft.yaml` by **0.158
THz**. The labels are correct.

**2. The two compounds use different cells.** Bi₂Te₃ is the 5-atom rhombohedral primitive
cell (15 branches); Sb₂Te₃ is the 15-atom conventional cell (45 branches). Each is
internally consistent — DFT and NEP always share a cell — but do not compare branch counts
across compounds.

## What the comparison shows

| Bi₂Te₃ | DFT | NEP | difference |
|---|---:|---:|---:|
| acoustic RMS error | — | — | 0.186 THz |
| optical RMS error | — | — | 0.149 THz (5.4 %) |
| phonon bandwidth | 4.215 | 3.807 THz | **−9.7 %** |
| optical manifold width | 2.980 | 2.504 THz | **−16.0 %** |
| acoustic bandwidth | 1.518 | — | — |

Sb₂Te₃ behaves the same way: −11.6 % bandwidth, −12.4 % manifold width.

The NEP reproduces the acoustic branches that carry the heat, but **systematically softens
the optical manifold**. This is discussed in the manuscript, including its consequence for
the acoustic+optical ↔ optical Umklapp phase space.

## A note on the M point

A supercell samples the dynamical matrix *exactly* only at wavevectors commensurate with
it. With `DIM = 3 3 2`, both Γ = (0,0,0) and K = (1/3,1/3,0) are commensurate, but
M = (1/2,0,0) is **not** — it requires an even in-plane expansion. M is therefore the only
one of the three high-symmetry points that is Fourier-interpolated, in the DFT curve as
much as in the NEP one. Recomputing with a 4×4×3 cell (four is even) removes the effect at
M entirely. This is why agreement at M looks worse than at K.
