# NEP training — which file reproduces the paper

**Use `nep.txt`.** It is the **generation-100,000** checkpoint, and it is the potential
used for *every* molecular-dynamics result in the paper.

| file | md5 | what it is |
|---|---|---|
| `nep.txt` | `545adcd43d8d8e19928a9f75937dee52` | **the paper's potential** (copy of the gen-100,000 checkpoint below) |
| `nep_y2026_m04_d08_h20_m24_s14_generation100000.txt` | `545adcd43d8d8e19928a9f75937dee52` | the same file, under its original name |
| `nep_y2026_m04_d09_h03_m51_s23_generation200000.txt` | `6f99078ae7409b5ce5b77c6581a0ed7c` | later checkpoint, **not used** |
| `nep_y2026_m04_d09_h11_m18_s39_generation300000.txt` | `c4e85ca2cef420247a9b89c9c28b0927` | end of training, **not used** |

All four files are 186,193 bytes, so file size does **not** distinguish them — check the md5.

> Earlier revisions of this repository shipped `nep.txt` as the generation-300,000
> checkpoint. Anyone who used that file would not reproduce the published results. It has
> been replaced by the generation-100,000 potential; the 300,000 checkpoint is still here
> under its own name.

## Why generation 100,000 and not 300,000

`nep.in` requests 300,000 generations and `loss.out` runs to the end, but the loss had
converged well before that. The generation-100,000 checkpoint was adopted for production.
The accuracy quoted in the paper is that checkpoint's, and the parity files here reproduce
it exactly:

| | energy | force | 
|---|---|---|
| train (`energy_train.out`, `force_train.out`) | 2.25 meV/atom | 103.44 meV/Å |
| test (`energy_test.out`, `force_test.out`) | 2.40 meV/atom | 105.88 meV/Å |

## Training settings

`nep.in` is the complete, unmodified input. The loss weights are:

    lambda_e      1.0
    lambda_f      1.0
    lambda_v      0.8
    lambda_shear  5.0

`lambda_shear` is a real NEP keyword (GPUMD `src/main_nep/parameters.cu`); the shear
weighting used here is the value above.

Note these are **per-property** weights inside the loss. They are distinct from the
**per-structure** weight of 1.0 applied uniformly to every structure in the training
database — the two are independent.

## Other files

| file | what |
|---|---|
| `nep.in` | training input, as run |
| `loss.out` | loss history to generation 300,000 |
| `energy_*.out`, `force_*.out` | parity data for the gen-100,000 checkpoint |
| `*.restart` | restart files for the corresponding checkpoints |
| `fig/` | loss and parity figures |
