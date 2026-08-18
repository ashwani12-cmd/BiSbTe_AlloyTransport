from ase.io import read
from ase.build import make_supercell
import numpy as np
import random

# ============================================================
# INPUT FILE
# ============================================================
input_file = "/home/IITB/multiscale-mechanics/ashwani12/ashwani/Bi2Te3_Sb2Te3/nep/elastic_const/energy_strain/dft/Bi2Te3/more_relax/scf/espresso.pwi"

# ============================================================
# SUPERCELL (IMPORTANT)
# Your unit cell = 15 atoms → this gives 60 atoms total
# ============================================================
SUPERCELL = (2, 2, 1)

# ============================================================
# TARGET COMPOSITIONS (best possible matches)
# Based on 24 Bi atoms
# ============================================================
cases = {
    "20pct": 5/24,   # 20.8%
    "40pct": 10/24,  # 41.7%
    "60pct": 14/24,  # 58.3%
    "80pct": 19/24   # 79.2%
}

# ============================================================
# FIXED SETTINGS
# ============================================================
SPECIES_ORDER = ["Bi", "Sb", "Te"]
MASSES = {"Bi": 208.980, "Sb": 121.760, "Te": 127.600}

random.seed(42)  # reproducibility

# ============================================================
# READ STRUCTURE
# ============================================================
atoms_primitive = read(input_file, format="espresso-in")
print(f"Primitive cell atoms: {len(atoms_primitive)}")

# ============================================================
# BUILD SUPERCELL
# ============================================================
nx, ny, nz = SUPERCELL
P = np.diag([nx, ny, nz])
atoms_base = make_supercell(atoms_primitive, P)

print(f"Supercell {nx}x{ny}x{nz}: {len(atoms_base)} atoms")

# ============================================================
# LOOP OVER COMPOSITIONS
# ============================================================
for tag, frac in cases.items():

    atoms = atoms_base.copy()
    symbols = atoms.get_chemical_symbols()

    # ---- Find Bi atoms ----
    bi_indices = [i for i, s in enumerate(symbols) if s == "Bi"]
    n_bi = len(bi_indices)

    n_sub = int(round(frac * n_bi))

    print(f"\n--- {tag} ---")
    print(f"Total Bi: {n_bi}, Sb replacing: {n_sub}")

    # ---- Random substitution ----
    replace_ids = random.sample(bi_indices, n_sub)

    for idx in replace_ids:
        symbols[idx] = "Sb"

    atoms.set_chemical_symbols(symbols)

    # ========================================================
    # BUILD ATOM TYPES
    # ========================================================
    type_map = {sp: i+1 for i, sp in enumerate(SPECIES_ORDER)}
    atom_types = [type_map[s] for s in symbols]

    # ========================================================
    # CELL → LAMMPS BOX
    # ========================================================
    cell = atoms.get_cell()
    pos  = atoms.get_positions()
    N    = len(atoms)

    a, b, c = cell[0], cell[1], cell[2]

    xhi = np.linalg.norm(a)
    xlo = 0.0
    xy  = np.dot(b, a / xhi)
    yhi = np.sqrt(np.dot(b, b) - xy**2)
    ylo = 0.0
    xz  = np.dot(c, a / xhi)
    yz  = (np.dot(b, c) - xy * xz) / yhi
    zhi = np.sqrt(max(np.dot(c, c) - xz**2 - yz**2, 0.0))
    zlo = 0.0

    is_triclinic = not (abs(xy) < 1e-8 and abs(xz) < 1e-8 and abs(yz) < 1e-8)

    # ---- Transform positions ----
    M = np.array([[xhi, 0,   0],
                  [xy,  yhi, 0],
                  [xz,  yz,  zhi]])

    frac_coords = np.linalg.solve(cell.T, pos.T).T
    lmp_pos = frac_coords @ M

    # ========================================================
    # OUTPUT FILE
    # ========================================================
    base = input_file.split("/")[-1].replace(".pwi", "")
    output_file = f"{base}_{nx}x{ny}x{nz}_{tag}.lmp"

    with open(output_file, "w") as f:
        f.write(f"# Bi2Te3 -> Sb alloy ({tag}) | {nx}x{ny}x{nz}\n\n")
        f.write(f"{N} atoms\n")
        f.write(f"3 atom types\n\n")

        f.write(f"{xlo:.6f} {xhi:.6f} xlo xhi\n")
        f.write(f"{ylo:.6f} {yhi:.6f} ylo yhi\n")
        f.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n")
        if is_triclinic:
            f.write(f"{xy:.6f} {xz:.6f} {yz:.6f} xy xz yz\n")
        f.write("\n")

        f.write("Masses\n\n")
        for i, sp in enumerate(SPECIES_ORDER):
            f.write(f"{i+1} {MASSES[sp]:.3f}  # {sp}\n")
        f.write("\n")

        f.write("Atoms  # atomic\n\n")
        for i in range(N):
            t = atom_types[i]
            x, y, z = lmp_pos[i]
            f.write(f"{i+1} {t} {x:.6f} {y:.6f} {z:.6f}\n")

    print(f"Written: {output_file}")

print("\n✅ DONE: All alloy structures generated correctly!")
