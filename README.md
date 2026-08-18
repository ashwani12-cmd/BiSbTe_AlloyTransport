# BiSbTe_AlloyTransport

Thermal transport and lattice dynamics of Bi–Sb–Te alloys using neural equivariant potentials (NEP) and molecular dynamics simulations.

This repository contains the DFT reference structures, training dataset, trained NEP model, and analysis scripts used to study lattice thermal transport and mechanical properties in (Bi₁₋ₓSbₓ)₂Te₃ alloys and Bi₂Te₃/Sb₂Te₃ interfaces — including bulk thermal conductivity (κ), interfacial Kapitza resistance (R_K), and elastic constants.

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
│   └── superlattice_NEMD_singleIF_FIXED/     # Single-interface Kapitza resistance NEMD runs
│       ├── SI_10Bi_10Sb_TeTe_10x14/          # 10+10 QL, 10x14 in-plane, Te-Te termination
│       ├── SI_15Bi_15Sb_TeTe_10x14/          # 15+15 QL
│       ├── SI_30Bi_30Sb_TeTe_10x14/          # 30+30 QL
│       ├── analyze_interface.py              # Extracts R_K, G_K, temperature profiles from compute.out
│       ├── plot_rk_three_fig.py              # Plots R_K vs L, dT_interface vs L, finite-size correction
│       ├── Fig1_RK_vs_L.pdf / .png
│       ├── Fig2_dT_vs_L.pdf / .png
│       ├── Fig3_FSC.pdf / .png
│       ├── RK_summary_JPz.pdf / .png
│       └── T_profile_paper_SI_*.pdf / .png   # Steady-state temperature profiles per configuration
│
└── elastic_constants_LAMMPS/
    ├── pwi_lmp_alloy.py              # Converts QE-relaxed alloy structures to LAMMPS data files
    ├── espresso_Bi2Te3.pwi           # DFT-relaxed Bi2Te3 structure (input to the converter)
    ├── Bi2Te3.lmp                    # Bi2Te3 LAMMPS data file
    ├── Bi2Te3/                       # Elastic constant calculation, pure Bi2Te3
    ├── Sb2Te3/                       # Elastic constant calculation, pure Sb2Te3
    │   ├── Sb2Te3.lmp
    │   ├── in.elastic                # LAMMPS elastic-constant driver script
    │   ├── init.mod, potential.mod, displace.mod   # LAMMPS include scripts (setup, potential, strain)
    │   ├── log.lammps, out.dat, job.*.out/.err      # Run log and results
    │   └── submit.sh
    └── BiSbTe{20,40,60,80}/          # Elastic constants for each alloy composition (20-80% Sb)
        ├── espresso_2x2x1_*pct.lmp   # 2x2x1 supercell LAMMPS data file for this composition
        ├── in.elastic, init.mod, potential.mod, displace.mod
        ├── restart.equil             # Equilibrated restart configuration
        ├── nep_y2026_*.txt           # NEP potential file used for this run
        ├── log.lammps, out.dat, job.*.out/.err
        └── submit.sh
```

## Workflow overview

1. **Structure generation** (`BiSbTe_alloy_endpoints/`) — DFT-relaxed primitive and conventional cells for Bi₂Te₃, Sb₂Te₃, and (Bi₁₋ₓSbₓ)₂Te₃ alloy compositions (20/40/60/80% Sb) computed in Quantum ESPRESSO with LDA ONCV pseudopotentials.
2. **Training set construction** (`training_data_json_format/`) — Equilibrium structures are perturbed via random displacement, uniaxial strain, and shear strain to sample the configuration space needed for a robust interatomic potential, without requiring full AIMD. The dataset is converted to NEP and MTP formats (`training_data_nep_format/`, `training_data_mtp_format/`) via `script/json_to_nepxyz.ipynb`.
3. **NEP training** (`nep_train/`) — A neural equivariant potential is trained on energies, forces, virials, and stresses from the dataset, with parity plots (`fig/`) and the loss curve used to validate accuracy on held-out test data.
4. **Bulk thermal conductivity** (`bulk_kappa_NEMD/`) — NEMD simulations of pure Bi₂Te₃ along the X, Y, and Z crystallographic directions, run over a series of system lengths for finite-size extrapolation of κ.
5. **Interfacial thermal transport** (`interface_kappa_NEMD/`) — Single-interface Bi₂Te₃/Sb₂Te₃ NEMD simulations extract the Kapitza resistance (R_K) and interfacial thermal conductance (G_K) across multiple system lengths, serving as an out-of-distribution transferability test of the trained NEP (no interface configurations were included in training).
6. **Elastic constants** (`elastic_constants_LAMMPS/`) — QE-relaxed structures for Bi₂Te₃, Sb₂Te₃, and each alloy composition are converted to LAMMPS data files (`pwi_lmp_alloy.py`) and run through LAMMPS' standard `in.elastic` strain-displacement workflow (using the trained NEP as the interatomic potential) to extract the full elastic constant tensor for each composition.

## Requirements

- [GPUMD](https://github.com/brucefan1983/GPUMD) for NEMD simulations and NEP training
- [LAMMPS](https://www.lammps.org/) (with NEP pair-style support) for elastic constant calculations
- [Quantum ESPRESSO](https://www.quantum-espresso.org/) for DFT reference calculations
- Python 3 with `numpy`, `matplotlib`, and `ase` for structure generation and post-processing

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

## Citation

Ashwani Kushwaha, Abhishek Kumar, Amit Singh. "Atomistic Study of Alloying Effects in (Bi₁₋ₓSbₓ)₂Te₃ Thermoelectric Materials Using Neural Equivariant Potentials." Department of Mechanical Engineering, IIT Bombay, Mumbai 400076, India. (Manuscript in preparation.)
