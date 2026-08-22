#!/usr/bin/env python3
"""
Generate QE strain inputs for all 7 distortion sets (hexagonal crystal).
Output structure:
  C11/
    delta_-0.004/espresso.pwi
    delta_-0.002/espresso.pwi
    delta_+0.000/espresso.pwi
    delta_+0.002/espresso.pwi
    delta_+0.004/espresso.pwi
  C12/
    ...
  C33/ C44/ C13p/ C13m/ Bulk/
"""

import os
import numpy as np
from ase.io import read, write
from ase.data import atomic_masses

# ============================================================
# USER SETTINGS
# ============================================================
input_file   = "/home/IITB/multiscale-mechanics/ashwani12/ashwani/Bi2Te3_Sb2Te3/nep/elastic_const/energy_strain/dft/Bi2Te3/more_relax/scf/espresso.pwi"
strain_range = np.arange(-0.004, 0.0041, 0.002)
pseudo_dir   = '/home/IITB/multiscale-mechanics/amit.k.singh/PP/PBE_ONCV/'
pseudos      = {"Bi": "Bi.UPF", "Te": "Te.UPF"}
kpts         = (12, 12, 1)

# ============================================================
# 7 STRAIN SETS
# ============================================================
STRAIN_SETS = {
    "C11" : lambda d: np.array([[1+d, 0,   0  ],
                                 [0,   1+d, 0  ],
                                 [0,   0,   1  ]]),

    "C12" : lambda d: np.array([[1+d, 0,   0  ],
                                 [0,   1-d, 0  ],
                                 [0,   0,   1  ]]),

    "C33" : lambda d: np.array([[1,   0,   0  ],
                                 [0,   1,   0  ],
                                 [0,   0,   1+d]]),

    "C44" : lambda d: np.array([[1,   0,   0  ],
                                 [0,   1,   d  ],
                                 [0,   d,   1  ]]),

    "C13p": lambda d: np.array([[1+d, 0,   0  ],
                                 [0,   1,   0  ],
                                 [0,   0,   1+d]]),

    "C13m": lambda d: np.array([[1+d, 0,   0  ],
                                 [0,   1,   0  ],
                                 [0,   0,   1-d]]),

    "Bulk": lambda d: np.array([[1+d, 0,   0  ],
                                 [0,   1+d, 0  ],
                                 [0,   0,   1+d]]),
}

# ============================================================
# QE INPUT BUILDER
# ============================================================
def prepare_input_data(atoms):
    elements = sorted(set(atoms.get_chemical_symbols()))
    atomic_species = []
    for el in elements:
        idx = atoms.get_chemical_symbols().index(el)
        mass = atomic_masses[atoms[idx].number]
        atomic_species.append([el, mass, pseudos[el]])
    return {
        "control": {
            "calculation"  : "relax",
            "verbosity"    : "high",
            "restart_mode" : "from_scratch",
            "etot_conv_thr": 1e-10,
            "forc_conv_thr": 1e-10,
            "tstress"      : True,
            "tprnfor"      : True,
            "outdir"       : "./",
            "disk_io"      : "none",
            "pseudo_dir"   : pseudo_dir,
        },
        "system": {
            "ibrav"          : 0,
            "ecutwfc"        : 120.0,
            "occupations"    : "smearing",
            "smearing"       : "mp",
            "degauss"        : 0.01,
            "nat"            : len(atoms),
            "ntyp"           : len(elements),
            "vdw_corr"       : "grimme-d3",
            "dftd3_version"  : 4,
            "dftd3_threebody": True,
        },
        "electrons": {
            "electron_maxstep": 400,
            "conv_thr"        : 1e-12,
            "mixing_beta"     : 0.3,
        },
        "atomic_species": atomic_species,
    }

# ============================================================
# READ STRUCTURE
# ============================================================
atoms = read(input_file, format="espresso-in")
cell0 = np.array(atoms.cell)
print(f"📘 Loaded {len(atoms)} atoms from {input_file}")
print("Initial cell (Å):\n", cell0)

# ============================================================
# LOOP: 7 folders × 5 strains
# ============================================================
for set_name, F_func in STRAIN_SETS.items():
    set_dir = set_name                          # e.g.  C11/
    os.makedirs(set_dir, exist_ok=True)
    print(f"\n📁 {set_dir}/")

    for d in strain_range:
        sub_dir = os.path.join(set_dir, "relax", f"delta_{d:+.3f}")   # e.g.  C11/relax/delta_-0.004/
        os.makedirs(sub_dir, exist_ok=True)

        F        = F_func(d)
        new_cell = cell0 @ F

        atoms_s  = atoms.copy()
        atoms_s.set_cell(new_cell, scale_atoms=True)

        input_data = prepare_input_data(atoms_s)

        pwi_path = os.path.join(sub_dir, "espresso.pwi")
        write(
            pwi_path,
            atoms_s,
            format           = "espresso-in",
            input_data       = input_data,
            pseudopotentials = {a[0]: a[2] for a in input_data["atomic_species"]},
            kpts             = kpts,
        )

        with open(pwi_path, "r+") as f:
            content = f.read()
            f.seek(0, 0)
            f.write(f"! Set={set_name}  delta={d:+.6f}\n" + content)

        # ---- write submit.sh ----
        submit_path = os.path.join(sub_dir, "submit.sh")
        with open(submit_path, "w") as f:
            f.write(f"""#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --partition=debug
#SBATCH --output={set_name}-%A_%03a.out
#SBATCH --error={set_name}-%A_%03a.err
#SBATCH --job-name={set_name}_{d:+.3f}
#SBATCH --mem=800

# Load environment
module load openmpi/4.1.4
module load nvhpc/24.11
module load intel-oneapi-mkl/2024.2.2-oneapi-2025.0.1-5u4sz3m

# Set OpenMP threads
export OMP_NUM_THREADS=1

# Path to correct mpirun (HPC-X inside nvhpc)
MPIRUN_PATH="/home/apps/hpc_sdk/Linux_x86_64/24.11/comm_libs/12.6/hpcx/hpcx-2.20/ompi/bin/mpirun"

# Path to your pw.x binary
PWX_BIN="/home/IITB/multiscale-mechanics/amit.k.singh/software/builds/qe-7.1/build/bin/pw.x"

# Move to job submission directory
cd $SLURM_SUBMIT_DIR

# Run Quantum ESPRESSO
time ${{MPIRUN_PATH}} -np "$SLURM_NTASKS" ${{PWX_BIN}} < espresso.pwi > espresso.pwo
""")

        print(f"    📄 {sub_dir}/espresso.pwi  +  submit.sh")

print(f"\n✅ Done — 7 folders, 5 strains each.")
