# BiSbTe_AlloyTransport

Thermal transport and lattice dynamics of Bi–Sb–Te alloys using neural equivariant potentials (NEP) and molecular dynamics simulations.

This repository contains the training data, trained potential, and analysis scripts used to study lattice thermal transport in (Bi₁₋ₓSbₓ)₂Te₃ alloys and Bi₂Te₃/Sb₂Te₃ interfaces, including bulk thermal conductivity (κ) and interfacial Kapitza resistance (R_K) via non-equilibrium molecular dynamics (NEMD).

## Repository structure

```
BiSbTe_AlloyTransport/
├── BiSbTe_alloy_endpoints/       # DFT-relaxed structures for the alloy end-members and compositions
├── nep_train/                    # NEP training run: inputs, checkpoints, and parity plots
│   ├── nep.in                    # NEP training configuration
│   ├── submit.sh                 # Cluster job submission script
│   ├── nep_y2026_*.txt/.restart  # Trained NEP model checkpoints (100k / 200k / 300k generations)
│   ├── energy_train.out, energy_test.out
│   ├── force_train.out,  force_test.out
│   ├── stress_train.out, stress_test.out
│   ├── virial_train.out, virial_test.out
│   ├── loss.out                  # Training loss history
│   └── fig/                      # Parity plots (energy/force/stress, train & test) and loss curve
├── training_data_nep_format/     # Training dataset in NEP (extxyz) format
├── training_data_mtp_format/     # Same dataset converted to MTP format
├── training_data_json_format/    # Same dataset converted to JSON format
├── bulk_kappa_NEMD/
│   └── Bi2Te3/                   # Bulk Bi2Te3 thermal conductivity via NEMD
│       ├── X/, Y/, Z/            # Transport direction (crystallographic axis)
│       │   ├── L_0250A/ .. L_1000A/   # System-length series for finite-size extrapolation
│       │   ├── JP_vs_bin_JP*.pdf
│       │   ├── finite_size_JP*.pdf
│       │   ├── temperature_profiles_JP*.pdf
│       │   └── post_processing_nemd_Jp.py
├── interface_kappa_NEMD/
│   ├── create_superlattices_FIXED.py       # Builds strained, stacked Bi2Te3/Sb2Te3 interface supercells
│   ├── espresso_Bi2Te3_conventional.pwi    # DFT-relaxed Bi2Te3 conventional cell (QE input)
│   ├── espresso_Sb2Te3_conventional.pwi    # DFT-relaxed Sb2Te3 conventional cell (QE input)
│   └── superlattice_NEMD_singleIF_FIXED/   # Single-interface Kapitza resistance NEMD runs
│       ├── SI_10Bi_10Sb_TeTe_10x14/        # 10+10 QL, 10x14 in-plane, Te-Te termination
│       ├── SI_15Bi_15Sb_TeTe_10x14/        # 15+15 QL
│       ├── SI_30Bi_30Sb_TeTe_10x14/        # 30+30 QL
│       ├── analyze_interface.py            # Extracts R_K, G_K, temperature profiles from compute.out
│       ├── plot_rk_three_fig.py            # Plots R_K vs L, dT_interface vs L, finite-size correction
│       ├── Fig1_RK_vs_L.pdf/.png
│       ├── Fig2_dT_vs_L.pdf/.png
│       ├── Fig3_FSC.pdf/.png
│       ├── RK_summary_JPz.pdf/.png
│       └── T_profile_paper_SI_*.pdf/.png   # Steady-state temperature profiles per configuration
└── script/                       # Miscellaneous helper/utility scripts
```

## Workflow overview

1. **Structure generation** — DFT-relaxed conventional cells (Quantum ESPRESSO) for Bi₂Te₃, Sb₂Te₃, and (Bi₁₋ₓSbₓ)₂Te₃ alloy compositions are used as the basis for supercell construction. Interface supercells are built by averaging the in-plane lattice constants of the two end-members and stacking whole quintuple layers (QLs) along the transport direction, with fixed and thermostatted regions arranged to bury the periodic-image interface.
2. **NEP training** (`nep_train/`) — A neural equivariant potential is trained on energies, forces, virials, and stresses from the DFT dataset (`training_data_*_format/`), with parity plots and loss curves used to validate accuracy on held-out test data.
3. **Bulk thermal conductivity** (`bulk_kappa_NEMD/`) — NEMD simulations of pure Bi₂Te₃ along the X, Y, and Z crystallographic directions, run over a series of system lengths for finite-size extrapolation of κ.
4. **Interfacial thermal transport** (`interface_kappa_NEMD/`) — Single-interface Bi₂Te₃/Sb₂Te₃ NEMD simulations to extract the Kapitza resistance (R_K) and interfacial thermal conductance (G_K), evaluated across multiple system lengths as an out-of-distribution transferability test of the trained NEP (no interface configurations were included in training).

## Requirements

- [GPUMD](https://github.com/brucefan1983/GPUMD) for NEMD simulations and NEP training
- [Quantum ESPRESSO](https://www.quantum-espresso.org/) for DFT reference calculations
- Python 3 with `numpy`, `matplotlib`, and `ase` for structure generation and post-processing

## Citation

If you use this code or data, please cite the associated manuscript (in preparation).
