from ase.io import read, write
from ase.build import make_supercell


# =========================
# INPUT FILE (CHANGE PATH)
# =========================

input_file = "/home/IITB/multiscale-mechanics/ashwani12/ashwani/Bi2Te3_Sb2Te3/conventional_cell/Bi2Te3/vc-relax/espresso.pwo"


# =========================
# READ STRUCTURE
# =========================

atoms = read(input_file, format="espresso-out")


# =========================
# BUILD 6x6x6 SUPERCELL
# =========================

supercell_matrix = [
    [6, 0, 0],
    [0, 6, 0],
    [0, 0, 6]
]

sc = make_supercell(atoms, supercell_matrix)


# =========================
# WRITE TEMP LAMMPS FILE
# =========================

write("lammps_tmp.data", sc, format="lammps-data", atom_style="atomic")


# =========================
# FIX TYPES: FORCE 3 TYPES
# =========================

infile = "lammps_tmp.data"
outfile = "lammps.data"

with open(infile) as f:
    lines = f.readlines()


new_lines = []
skip_masses = False

for line in lines:

    # Fix header
    if "atom types" in line:
        new_lines.append("3 atom types\n")
        continue


    # Rewrite Masses
    if line.strip() == "Masses":
        new_lines.append("Masses\n\n")
        new_lines.append("1 208.98   # Bi\n")
        new_lines.append("2 121.76   # Sb\n")
        new_lines.append("3 127.60   # Te\n\n")
        skip_masses = True
        continue


    # Skip old Masses
    if skip_masses:
        if line.strip() == "":
            skip_masses = False
        continue


    new_lines.append(line)


final_lines = []

for line in new_lines:

    parts = line.split()

    # Fix atom types in Atoms section
    if len(parts) >= 5 and parts[0].isdigit():

        old_type = parts[1]

        # Old: 1=Bi , 2=Te
        if old_type == "1":
            parts[1] = "1"   # Bi
        elif old_type == "2":
            parts[1] = "3"   # Te

        line = "  ".join(parts) + "\n"

    final_lines.append(line)


with open(outfile, "w") as f:
    f.writelines(final_lines)


print("================================")
print("DONE ✅")
print("Output file:", outfile)
print("Type mapping: 1=Bi, 2=Sb(dummy), 3=Te")
print("================================")

