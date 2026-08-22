# BiSbTe_AlloyTransport

Thermal transport and lattice dynamics of Bi–Sb–Te alloys using neural equivariant potentials (NEP) and molecular dynamics simulations.

This repository contains the DFT reference structures, training dataset, trained NEP model, and analysis scripts used to study lattice thermal transport and mechanical properties in (Bi₁₋ₓSbₓ)₂Te₃ alloys and Bi₂Te₃/Sb₂Te₃ interfaces — including bulk thermal conductivity (κ), interfacial Kapitza resistance (R_K), and elastic constants (Cij).

## Repository structure

```
BiSbTe_AlloyTransport/
├── BiSbTe_alloy_endpoints/
│   ├── Bi.UPF, Sb.UPF, Te.UPF                        # LDA ONCV pseudopotentials (Quantum ESPRESSO)
│   ├── espresso_Bi2Te3_primitive.pwi                 # DFT-relaxed primitive cells, end-members
│   ├── espresso_Bi2Te3_conventional.pwi
│   ├── espresso_Bi2Te3_332_supercell_primitive.pwi
│   ├── espresso_Sb2Te3_primitive.pwi
│   ├── espresso_Sb2Te3_conventional.pwi
│   └── espresso_BiSbTe{20,40,60,80}_primitive.pwi    # Alloy compositions, 20-80% Sb substitution
│
├── training_data_json_format/
│   ├── README.md                    # Detailed per-file naming convention and sampling counts
│   ├── equib_structure/             # DFT-relaxed equilibrium structures per composition
│   ├── near_harmonic/               # Near-equilibrium / small-perturbation configurations
│   ├── random_disp_0.3A/            # Random atomic displacements up to 0.3 A (200 configs/composition)
│   ├── shear_strain/                # 1% and 3% shear, with/without 0.3 A perturbation
│   └── uniaxial_strain/             # Uniaxial strain along X/Y/Z, -10% to +10% (11 points each)
├── training_data_nep_format/        # Same dataset converted to NEP/extxyz (train.xyz, test.xyz)
├── training_data_mtp_format/        # Same dataset converted to MTP format (train.cfg, test.cfg)
├── script/
│   └── json_to_nepxyz.ipynb         # Converts training_data_json_format -> NEP extxyz format
│
├── nep_train/
│   ├── nep.in                       # NEP training configuration
│   ├── submit.sh                    # Cluster job submission script
│   ├── train.xyz, test.xyz          # Training/test sets used for this run
│   ├── nep_y2026_*.txt / .restart   # Trained NEP checkpoints (100k / 200k / 300k generations)
│   ├── nep.txt, nep.restart         # Final/active model (symlink or copy of latest checkpoint)
│   ├── energy_train.out,  energy_test.out
│   ├── force_train.out,   force_test.out
│   ├── stress_train.out,  stress_test.out
│   ├── virial_train.out,  virial_test.out
│   ├── loss.out                     # Training loss history
│   ├── plot_loss.py                 # Plots loss curve from loss.out
│   ├── plot_parity_train.py         # Parity plots (train set): energy/force/stress
│   ├── plot_parity_test.py          # Parity plots (test set): energy/force/stress
│   └── fig/                         # Generated parity plots and loss curve (PDF + PNG, 100k)
│
├── bulk_kappa_NEMD/
│   └── Bi2Te3/                      # Bulk Bi2Te3 thermal conductivity via NEMD
│       ├── X/                       # Transport along the a-axis
│       ├── Y/                       # Transport along the b-axis
│       └── Z/                       # Transport along the c-axis (cross-plane)
│           ├── L_0250A/ .. L_1000A/      # System-length series for finite-size extrapolation
│           ├── JP_vs_bin_JP*.pdf
│           ├── finite_size_JP*.pdf
│           ├── temperature_profiles_JP*.pdf
│           └── post_processing_nemd_Jp.py
│
├── interface_kappa_NEMD/
│   ├── create_superlattices_FIXED.py         # Builds strained, stacked Bi2Te3/Sb2Te3 interface supercells
│   ├── espresso_Bi2Te3_conventional.pwi      # DFT-relaxed Bi2Te3 conventional cell (input to the builder)
│   ├── espresso_Sb2Te3_conventional.pwi      # DFT-relaxed Sb2Te3 conventional cell (input to the builder)
│   └── 200K/                                 # Single-interface NEMD at 200 K (T_avg=200 K, dT=30 K)
│   └── 300K/                                 # Single-interface NEMD at 300 K
│   └── 400K/                                 # Single-interface NEMD at 400 K
│   └── 500K/                                 # Single-interface NEMD at 500 K
│       ├── model.xyz                         # Input supercell: 30 QL Bi2Te3 + 30 QL Sb2Te3, 10x14 in-plane,
│       │                                     #   42,000 atoms, Lz ≈ 60–63 nm, A ≈ 22.3 nm²
│       │                                     #   Grouping method 0: 9 transport bins (G0=wall, G1=src,
│       │                                     #   G2–G7=transport, G8=snk) + method 1 (G1=Bi2Te3, G2=Sb2Te3)
│       ├── relaxed.xyz                       # Post-minimization + NPT-equilibrated structure
│       ├── run.in                            # GPUMD input: FIRE minimize → NPT (1 ns) → NEMD (5 ns)
│       │                                     #   heat_lan ensemble, Langevin thermostat on G1 (src) and G8 (snk)
│       │                                     #   compute_shc on grouping method 0, group 4 (last Bi2Te3 bin,
│       │                                     #   z = 221–300 Å): shc_out (ko) = spectral heat flux across
│       │                                     #   G4→G5 boundary = Bi2Te3|Sb2Te3 interface
│       ├── nemd_setup.txt                    # Human-readable geometry summary: group boundaries, Lz, A_cross,
│       │                                     #   interface location, reservoir bounds
│       ├── nep_y2026_m04_d08_h20_m24_s14_generation100000.txt   # NEP4 potential (Bi/Sb/Te, 100k generations)
│       ├── compute.out                       # Per-bin temperature + heat current (JP, JK) every 1 ps
│       │                                     #   shape: (5000, 65) — 8 T cols + 6×9 JP/JK cols
│       ├── compute_chunk.out                 # Fine-grained 1D temperature/density profile (104 bins, 1 ps output)
│       ├── thermo.out                        # Global thermodynamic output (T, KE, PE, stress, box) every 1 ps
│       ├── shc.out                           # Spectral heat current output from compute_shc:
│       │                                     #   499 rows: time-domain HCACF (ki, ko), t ∈ [−0.498, +0.498] ps
│       │                                     #   1000 rows: frequency-domain SHC, ω ∈ [0.05, 50] THz
│       ├── neighbor.out                      # Neighbor list statistics
│       ├── gpu_monitor.log                   # GPU utilisation log (A100, CUDA 12.8)
│       ├── gpumd.out / gpumd.err             # GPUMD stdout/stderr
│       ├── out.dat                           # GPUMD run summary (atom counts, group assignments, timing)
│       ├── analyze_interface.py              # Post-processing script: reads compute.out + relaxed.xyz,
│       │                                     #   fits linear T-profiles on each side of the interface,
│       │                                     #   extracts R_K, G_K, dT_interface, κ (global + per-side)
│       ├── T_profile_paper_{T}K.pdf/.png     # Publication-quality steady-state temperature profile
│       ├── JP_vs_bin_JPz.pdf/.png            # Per-bin heat flux JPz (all bins + fit bins highlighted)
│       └── JP_deviation_JPz.pdf/.png         # JPz deviation from fit-bin mean + block SEM error bars
│
│   Key results (single Bi2Te3|Sb2Te3 interface, Te-Te termination, 30+30 QL, 10×14 in-plane):
│
│   | T (K) | κ_Bi2Te3 (W/mK) | κ_Sb2Te3 (W/mK) | ΔT_IF (K) | R_K (×10⁻⁹ m²K/W) | G_K (MW/m²K) |
│   |-------|-----------------|-----------------|-----------|-------------------|--------------|
│   |  200  |     0.119       |     0.560       |   5.60    |       28.78       |    34.75     |
│   |  300  |     0.130       |     0.324       |   3.37    |       17.56       |    56.94     |
│   |  400  |     0.109       |     0.359       |   3.27    |       16.46       |    60.76     |
│   |  500  |     0.091       |     0.525       |   1.89    |       10.14       |    98.60     |
│
│   R_K decreases 65% from 200→500 K; G_K scales approximately as T^1.04,
│   consistent with diffuse mismatch model behaviour in the high-anharmonicity regime.
│
│   SHC spectral decomposition (shc_out, ko — heat flux across Bi2Te3|Sb2Te3 interface):
│   | T (K) | f_50% (THz) | f_80% (THz) | <2 THz | 2–5 THz | 5–10 THz | >10 THz |
│   |-------|-------------|-------------|--------|---------|----------|---------|
│   |  200  |    6.45     |    13.10    |  16.7% |  23.3%  |  28.7%   |  30.4%  |
│   |  300  |    6.60     |    13.70    |  16.3% |  22.7%  |  28.0%   |  31.9%  |
│   |  400  |    5.75     |    11.50    |  18.6% |  25.5%  |  29.7%   |  25.2%  |
│   |  500  |    5.35     |     9.90    |  19.5% |  27.2%  |  33.0%   |  19.2%  |
│
│   The median heat-carrying frequency (f_50%) redshifts from 6.6→5.4 THz with
│   increasing T, and the >10 THz optical contribution drops from 32%→19%,
│   consistent with anharmonic suppression of high-frequency modes at elevated T.
│
├── elastic_constants_LAMMPS/
│   ├── pwi_lmp_alloy.py              # Converts QE-relaxed alloy structures to LAMMPS data files
│   ├── espresso_Bi2Te3.pwi           # DFT-relaxed Bi2Te3 structure (input to the converter)
│   ├── Bi2Te3.lmp                    # Bi2Te3 LAMMPS data file
│   ├── Bi2Te3/                       # Elastic constant calculation (NEP potential), pure Bi2Te3
│   ├── Sb2Te3/                       # Elastic constant calculation (NEP potential), pure Sb2Te3
│   │   ├── Sb2Te3.lmp
│   │   ├── in.elastic                # LAMMPS elastic-constant driver script
│   │   ├── init.mod, potential.mod, displace.mod   # LAMMPS include scripts (setup, potential, strain)
│   │   ├── log.lammps, out.dat, job.*.out/.err      # Run log and results
│   │   └── submit.sh
│   └── BiSbTe{20,40,60,80}/          # Elastic constants for each alloy composition (20-80% Sb)
│       ├── espresso_2x2x1_*pct.lmp   # 2x2x1 supercell LAMMPS data file for this composition
│       ├── in.elastic, init.mod, potential.mod, displace.mod
│       ├── restart.equil             # Equilibrated restart configuration
│       ├── nep_y2026_*.txt           # NEP potential file used for this run
│       ├── log.lammps, out.dat, job.*.out/.err
│       └── submit.sh
│
└── elastic_constants_DFT_energy_strain/
    └── Bi2Te3/
        ├── generate_strain_all.py    # Generates the strained structures for every Cij set below
        ├── vc-relax/                 # Variable-cell relaxation to the equilibrium structure
        ├── scf/                      # Reference SCF at the unstrained (delta=0) equilibrium cell
        ├── Bulk/relax/delta_*/       # Isotropic-volume strain series (not currently used in the fit)
        ├── C11/                      # Strain set isolating C11+C12
        │   ├── relax/delta_*/            # Ionic relaxation at each strain delta
        │   └── scf/delta_*/              # Single-point SCF at each relaxed, strained delta
        ├── C12/                      # Strain set isolating C11-C12
        ├── C13m/                     # Strain set isolating (C11-2*C13+C33)/2
        ├── C13p/                     # Strain set isolating (C11+2*C13+C33)/2 (relax only, not yet used)
        ├── C33/                      # Strain set isolating C33/2
        ├── C44/                      # Strain set isolating 2*C44
        ├── get_Cij.py                # Solves the full Cij tensor from all 5 completed strain sets
        └── make_fit_figures.py       # Generates one labeled energy-vs-strain fit figure per set
```

## Workflow overview

1. **Structure generation** (`BiSbTe_alloy_endpoints/`) — DFT-relaxed primitive and conventional cells for Bi₂Te₃, Sb₂Te₃, and (Bi₁₋ₓSbₓ)₂Te₃ alloy compositions (20/40/60/80% Sb) computed in Quantum ESPRESSO with LDA ONCV pseudopotentials.
2. **Training set construction** (`training_data_json_format/`) — Equilibrium structures are perturbed via random displacement, uniaxial strain, and shear strain to sample the configuration space needed for a robust interatomic potential, without requiring full AIMD. The dataset is converted to NEP and MTP formats (`training_data_nep_format/`, `training_data_mtp_format/`) via `script/json_to_nepxyz.ipynb`.
3. **NEP training** (`nep_train/`) — A neural equivariant potential is trained on energies, forces, virials, and stresses from the dataset, with parity plots (`fig/`) and the loss curve used to validate accuracy on held-out test data.
4. **Bulk thermal conductivity** (`bulk_kappa_NEMD/`) — NEMD simulations of pure Bi₂Te₃ along the X, Y, and Z crystallographic directions, run over a series of system lengths for finite-size extrapolation of κ.
5. **Interfacial thermal transport** (`interface_kappa_NEMD/`) — Single-interface Bi₂Te₃/Sb₂Te₃ NEMD simulations at 200, 300, 400, and 500 K extract the Kapitza resistance (R_K) and interfacial thermal conductance (G_K). No interface configurations were included in training — these runs serve as a strict out-of-distribution transferability test of the NEP. Each temperature folder contains the full GPUMD run (FIRE minimisation → 1 ns NPT equilibration → 5 ns NEMD production on an A100 GPU), the spectral heat current output (`shc.out`), and publication-quality figures generated by `analyze_interface.py`.
6. **Elastic constants — NEP/LAMMPS** (`elastic_constants_LAMMPS/`) — QE-relaxed structures for Bi₂Te₃, Sb₂Te₃, and each alloy composition are converted to LAMMPS data files (`pwi_lmp_alloy.py`) and run through LAMMPS' standard `in.elastic` strain-displacement workflow (using the trained NEP as the interatomic potential) to extract the full elastic constant tensor for each composition.
7. **Elastic constants — DFT energy-strain** (`elastic_constants_DFT_energy_strain/`) — Independent DFT reference values for Bi₂Te₃, computed by applying small (±0.2%, ±0.4%) strains along five independent deformation modes, relaxing ions at fixed strained cell shape, and fitting the resulting `(E-E₀)/V₀` vs. strain `δ` curve to a quadratic. Each strain set's fit coefficient corresponds to a known linear combination of the Cij (see table below); solving the resulting 5×5 linear system yields the full independent set C11, C12, C13, C33, C44 (C66 follows algebraically as `(C11-C12)/2` for this trigonal symmetry, point group -3m — no separate strain set is needed for it).

### DFT elastic constants — Bi₂Te₃ results

Computed with `elastic_constants_DFT_energy_strain/Bi2Te3/get_Cij.py`:

| Strain set | Combination fit | Fit value (GPa) |
|---|---|---|
| C11  | C11 + C12                 | 92.52 |
| C12  | C11 − C12                 | 56.56 |
| C13m | (C11 − 2·C13 + C33) / 2   | 34.90 |
| C33  | C33 / 2                   | 24.01 |
| C44  | 2·C44                     | 68.72 |

Solved elastic constants:

| Constant | Value (GPa) |
|---|---|
| C11 | 74.54 |
| C12 | 17.98 |
| C13 | 26.38 |
| C33 | 48.02 |
| C44 | 34.36 |
| C66 | 28.28 *(= (C11−C12)/2, not independently fit)* |

`C11 > C33` reflects the expected anisotropy of a layered van der Waals material — Bi₂Te₃ is more compressible along the c-axis (cross-plane, van der Waals-bonded direction) than in-plane. Per-set fit quality and labeled energy-vs-strain curves (with the fit equation, coefficient, corresponding Cij combination, and R²) are generated by `make_fit_figures.py` and saved under `fit_figures/`.

**Usage:**
```bash
cd elastic_constants_DFT_energy_strain/Bi2Te3
python3 get_Cij.py .                          # prints all 6 Cij to stdout
python3 make_fit_figures.py . fit_figures     # saves one PDF+PNG fit figure per strain set
```

Note: the `C13p` and `Bulk` strain sets currently only have `relax/` inputs (no completed `scf/` outputs), so they are not used in the current fit — the 5 completed sets above are sufficient to solve the full Cij tensor without them.

## Requirements

- [GPUMD](https://github.com/brucefan1983/GPUMD) for NEMD simulations and NEP training
- [LAMMPS](https://www.lammps.org/) (with NEP pair-style support) for elastic constant calculations
- [Quantum ESPRESSO](https://www.quantum-espresso.org/) for DFT reference calculations
- Python 3 with `numpy`, `scipy`, `matplotlib`, and `ase` for structure generation and post-processing

## Notes for maintainers

- `training_data_json_format/README.md` documents the exact naming convention and configuration counts for each perturbation type (random displacement, shear, uniaxial strain) — see that file for details before adding new training data.
- A `del/` folder and a stray top-level `plot_parity_test.py` are currently still present in the repository from earlier scratch work and duplicate a file already tracked in `nep_train/`. These are scheduled for removal:

  ```bash
  git rm -r --cached del
  git rm --cached plot_parity_test.py
  echo "del/" >> .gitignore
  echo "plot_parity_test.py" >> .gitignore
  git add .gitignore
  git commit -m "Remove scratch del/ folder and duplicate plot_parity_test.py from tracking"
  git push
  ```
- `elastic_constants_DFT_energy_strain/` currently contains Bi₂Te₃ only. When Sb₂Te₃ and alloy-composition DFT elastic constants are added, nest them the same way (`elastic_constants_DFT_energy_strain/Sb2Te3/`, `.../BiSbTe40/`, etc.) to keep the structure consistent with `elastic_constants_LAMMPS/`.
- The `interface_kappa_NEMD/` temperature folders (200K–500K) each contain a full GPUMD run including `shc.out` for spectral heat current analysis. The `compute_shc` keyword in `run.in` uses `group 0 4` (grouping method 0, group ID 4 = last Bi₂Te₃ transport bin, z = 221–300 Å). The `shc_out (ko)` column in `shc.out` gives the spectral heat flux crossing the G4→G5 boundary, i.e. directly across the Bi₂Te₃|Sb₂Te₃ interface. Note: the 200K run predates the `compute_chunk` keyword addition and therefore has no `compute_chunk.out`; all other outputs are present and consistent.

## Citation

Ashwani Kushwaha, Abhishek Kumar, Amit Singh. "Atomistic Study of Alloying Effects in (Bi₁₋ₓSbₓ)₂Te₃ Thermoelectric Materials Using Neural Equivariant Potentials." Department of Mechanical Engineering, IIT Bombay, Mumbai 400076, India. (Manuscript in preparation.)
