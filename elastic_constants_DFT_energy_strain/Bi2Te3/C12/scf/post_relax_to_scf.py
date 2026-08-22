#!/usr/bin/env python3
"""
post_relax_to_scf.py
--------------------
Run this from:  C11/scf/
It reads relaxed coords from the hardcoded relax path and creates
delta_±X.XXX/ folders right here (in cwd) with espresso.pwi + submit.sh

Usage:
    cd .../test/C11/scf
    python post_relax_to_scf.py
"""

import os
import re

# ============================================================
# USER SETTINGS  — only change SET_NAME per strain set
# ============================================================
SET_NAME  = "C12"   # change to: C12, C33, C44, C13p, C13m, Bulk


RELAX_DIR = (
    "/home/IITB/multiscale-mechanics/ashwani12/ashwani/"
    "Bi2Te3_Sb2Te3/nep/elastic_const/energy_strain/dft/Bi2Te3/more_relax/"
    f"{SET_NAME}/relax"
)

SUBMIT_TEMPLATE = """\
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --partition=debug
#SBATCH --output={set_name}-%A_%03a.out
#SBATCH --error={set_name}-%A_%03a.err
#SBATCH --job-name={set_name}_{delta}_scf
#SBATCH --mem=800

module load openmpi/4.1.4
module load nvhpc/24.11
module load intel-oneapi-mkl/2024.2.2-oneapi-2025.0.1-5u4sz3m

export OMP_NUM_THREADS=1

MPIRUN_PATH="/home/apps/hpc_sdk/Linux_x86_64/24.11/comm_libs/12.6/hpcx/hpcx-2.20/ompi/bin/mpirun"
PWX_BIN="/home/IITB/multiscale-mechanics/amit.k.singh/software/builds/qe-7.1/build/bin/pw.x"

cd $SLURM_SUBMIT_DIR

time $MPIRUN_PATH -np $SLURM_NTASKS $PWX_BIN < espresso.pwi > espresso.pwo
"""

# ============================================================
# READ RELAXED POSITIONS FROM .pwo
# ============================================================
def read_relaxed_positions(pwo_path):
    with open(pwo_path, "r") as f:
        text = f.read()
    match = re.search(
        r"Begin final coordinates.*?ATOMIC_POSITIONS\s*\(angstrom\)\n(.*?)End final coordinates",
        text, re.DOTALL
    )
    if not match:
        raise RuntimeError(f"No 'Begin final coordinates' block in {pwo_path}")
    lines = match.group(1).strip().splitlines()
    return [re.sub(r"\s+", "  ", l.strip()) for l in lines if l.strip()]

# ============================================================
# BUILD SCF pwi
# ============================================================
def make_scf_pwi(relax_pwi, relaxed_positions, set_name, delta):
    with open(relax_pwi, "r") as f:
        lines = f.readlines()

    out = []
    skip_old_pos = False

    for line in lines:
        if line.startswith("! Set="):
            out.append(f"! Set={set_name}  delta={delta}  [SCF from relaxed]\n")
            continue

        if re.match(r"\s*calculation\s*=\s*'relax'", line):
            out.append(line.replace("'relax'", "'scf'"))
            continue

        if re.match(r"\s*ATOMIC_POSITIONS", line):
            out.append("ATOMIC_POSITIONS angstrom\n")
            for pos in relaxed_positions:
                out.append(pos + "\n")
            skip_old_pos = True
            continue

        if skip_old_pos:
            if re.match(r"\s*[A-Z_]{2,}", line) and not re.match(
                r"\s*(Bi|Te|Sb|Si|C |N |O |H |Mg)\s", line
            ):
                skip_old_pos = False
                out.append(line)
            continue

        out.append(line)

    return out

# ============================================================
# MAIN
# ============================================================
def main():
    scf_root = os.getcwd()   # folders created right here

    if not os.path.isdir(RELAX_DIR):
        print(f"ERROR: RELAX_DIR not found:\n   {RELAX_DIR}")
        return

    delta_folders = sorted(
        d for d in os.listdir(RELAX_DIR)
        if d.startswith("delta_") and os.path.isdir(os.path.join(RELAX_DIR, d))
    )

    if not delta_folders:
        print(f"ERROR: No delta_* folders in {RELAX_DIR}")
        return

    print(f"Reading relax : {RELAX_DIR}")
    print(f"Writing scf   : {scf_root}\n")

    for delta in delta_folders:
        relax_sub = os.path.join(RELAX_DIR, delta)
        pwo_path  = os.path.join(relax_sub, "espresso.pwo")
        pwi_path  = os.path.join(relax_sub, "espresso.pwi")

        if not os.path.isfile(pwo_path):
            print(f"  SKIP  {delta}  — espresso.pwo missing")
            continue
        if not os.path.isfile(pwi_path):
            print(f"  SKIP  {delta}  — espresso.pwi missing")
            continue

        try:
            relaxed_pos = read_relaxed_positions(pwo_path)
        except RuntimeError as e:
            print(f"  ERROR {delta}  — {e}")
            continue

        scf_lines = make_scf_pwi(pwi_path, relaxed_pos, SET_NAME, delta)

        out_dir = os.path.join(scf_root, delta)
        os.makedirs(out_dir, exist_ok=True)

        with open(os.path.join(out_dir, "espresso.pwi"), "w") as f:
            f.writelines(scf_lines)

        with open(os.path.join(out_dir, "submit.sh"), "w") as f:
            f.write(SUBMIT_TEMPLATE.format(set_name=SET_NAME, delta=delta))

        print(f"  OK  {delta}/   ({len(relaxed_pos)} atoms)")

    print(f"\nDone — {len(delta_folders)} SCF folders created.")

if __name__ == "__main__":
    main()
