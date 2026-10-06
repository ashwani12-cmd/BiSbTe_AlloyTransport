import os
import shutil
import numpy as np
from ase.io import read, write
from ase.build import make_supercell
from collections import Counter

# ============================================================
# USER SETTINGS
# ============================================================
BI2TE3_PWI  = "./espresso_Bi2Te3.pwi"
SB2TE3_PWI  = "./espresso_Sb2Te3.pwi"

import glob as _glob
_nep_files = _glob.glob("./nep_y*.txt") + _glob.glob("./nep*.txt")
if len(_nep_files) == 0:
    raise FileNotFoundError("No NEP potential file found!")
elif len(_nep_files) > 1:
    _nep_files.sort(reverse=True)
    print(f"  WARNING: Multiple NEP files found, using: {_nep_files[0]}")
POTENTIAL = _nep_files[0]
print(f"  Auto-detected NEP potential : {POTENTIAL}")

_submit_files = _glob.glob("./submit.sh") + _glob.glob("./*.sh")
SUBMIT_SH = _submit_files[0] if _submit_files else None
if SUBMIT_SH: print(f"  Auto-detected submit script : {SUBMIT_SH}")

TRANSPORT_DIRS = ['x', 'y', 'z']
TARGET_LENGTHS = [25.0, 50.0, 75.0, 100.0]   # nm
SB_FRACTIONS   = [0.2, 0.4, 0.6, 0.8]
RANDOM_SEED    = 1

T_TARGET        = 300
T_DELTA         = 20
THERMO_COUPLING = 100
MAX_OMEGA       = 50.0
EQ_STEPS        = 1000000
PROD_STEPS      = 5000000

TARGET_A_CROSS = 120.0   # nm² — FIXED cross-section target for ALL directions/materials
# MAX_ATOMS removed — no cap

EC_Bi2Te3 = dict(C11=55.0, C22=55.0, C33=10.0, C44=8.0,  C55=8.0,  C66=20.0)
EC_Sb2Te3 = dict(C11=63.0, C22=63.0, C33=14.0, C44=10.0, C55=10.0, C66=24.0)

FRAC_WALL = 0.02
FRAC_SRC  = 0.18
FRAC_SNK  = 0.18
N_MID     = 6
G_SRC     = 1
G_SNK     = 8

DIR_MAP = {'x': 0, 'y': 1, 'z': 2}

# ============================================================
# HELPERS
# ============================================================

def interp_ec(sb_frac):
    return {k: (1-sb_frac)*EC_Bi2Te3[k] + sb_frac*EC_Sb2Te3[k]
            for k in EC_Bi2Te3}


def get_uc_cross_area(atoms, transport_dir):
    """
    Compute the cross-sectional area of the unit cell
    for the two axes perpendicular to transport_dir.
    """
    cell = atoms.cell[:]
    a, b, c = cell[0], cell[1], cell[2]
    if   transport_dir == 'x': return np.linalg.norm(np.cross(b, c))
    elif transport_dir == 'y': return np.linalg.norm(np.cross(a, c))
    else:                      return np.linalg.norm(np.cross(a, b))


def build_supercell(atoms, transport_dir, target_L_nm, fixed_cross_reps=None):
    """
    Build supercell.
    - Transport reps scale with target_L_nm.
    - Cross reps are fixed to hit TARGET_A_CROSS (120 nm²).
    - fixed_cross_reps: pass previously computed dict to reuse; None = compute now.
    Returns: sc, nx, ny, nz, cross_reps
    """
    uc      = atoms.cell.lengths()
    dirs    = ['x', 'y', 'z']
    dir_idx = DIR_MAP[transport_dir]
    cross_dirs = [d for d in dirs if d != transport_dir]

    # Transport reps — changes per target_L
    uc_t     = uc[dir_idx]
    min_reps = 24 if uc_t < 10.0 else 1
    nx_t     = max(int(np.ceil((target_L_nm * 10) / uc_t)), min_reps)

    # Cross reps — solve for TARGET_A_CROSS
    if fixed_cross_reps is None:
        # Unit cell cross area (Å²) with 1×1 cross reps
        uc_cross_A2 = get_uc_cross_area(atoms, transport_dir)  # Å²
        target_A2   = TARGET_A_CROSS * 100.0                   # nm² → Å²

        # We need n_cross² × uc_cross_A2 ≈ target_A2
        # (assumes the two cross dims have equal reps, which is true for hex/rhombohedral)
        # More generally, solve n = ceil(sqrt(target/uc)) per cross axis
        # But since both cross dims are equal for Bi2Te3/Sb2Te3, n_a = n_b = n
        n_cross_float = np.sqrt(target_A2 / uc_cross_A2)
        n_cross = max(int(np.ceil(n_cross_float)), 3)

        cross_reps = {cd: n_cross for cd in cross_dirs}
    else:
        cross_reps = fixed_cross_reps

    if   transport_dir == 'x': nx, ny, nz = nx_t, cross_reps['y'], cross_reps['z']
    elif transport_dir == 'y': nx, ny, nz = cross_reps['x'], nx_t, cross_reps['z']
    else:                      nx, ny, nz = cross_reps['x'], cross_reps['y'], nx_t

    sc = make_supercell(atoms, [[nx,0,0],[0,ny,0],[0,0,nz]])
    return sc, nx, ny, nz, cross_reps


def assign_groups(supercell, transport_axis):
    L  = supercell.cell.lengths()[transport_axis]
    fp = supercell.get_scaled_positions()[:, transport_axis] * L

    frac_mid = (1.0 - 2*FRAC_WALL - FRAC_SRC - FRAC_SNK) / N_MID
    bf = [0.0, FRAC_WALL, FRAC_WALL + FRAC_SRC]
    for i in range(N_MID):
        bf.append(FRAC_WALL + FRAC_SRC + (i+1)*frac_mid)
    bf += [1.0 - FRAC_WALL, 1.0]
    zb = [b * L for b in bf]

    gids = np.zeros(len(supercell), dtype=int)
    for i, p in enumerate(fp):
        if   p < zb[1] or p >= zb[9]: gids[i] = 0
        elif p < zb[2]:                gids[i] = 1
        elif p < zb[3]:                gids[i] = 2
        elif p < zb[4]:                gids[i] = 3
        elif p < zb[5]:                gids[i] = 4
        elif p < zb[6]:                gids[i] = 5
        elif p < zb[7]:                gids[i] = 6
        elif p < zb[8]:                gids[i] = 7
        elif p < zb[9]:                gids[i] = 8
        else:                          gids[i] = 0
    return gids, zb, L


def get_z_bounds(L_transport):
    frac_mid = (1.0 - 2*FRAC_WALL - FRAC_SRC - FRAC_SNK) / N_MID
    bf = [0.0, FRAC_WALL, FRAC_WALL + FRAC_SRC]
    for i in range(N_MID):
        bf.append(FRAC_WALL + FRAC_SRC + (i+1)*frac_mid)
    bf += [1.0 - FRAC_WALL, 1.0]
    return [b * L_transport for b in bf]


def get_A_cross(supercell, transport_dir):
    cell = supercell.cell[:]
    a, b, c = cell[0], cell[1], cell[2]
    if   transport_dir == 'x': return np.linalg.norm(np.cross(b, c))
    elif transport_dir == 'y': return np.linalg.norm(np.cross(a, c))
    else:                      return np.linalg.norm(np.cross(a, b))


def substitute_bi_sb(supercell, sb_frac, seed):
    syms = np.array(supercell.get_chemical_symbols())
    bi_idx = np.where(syms == 'Bi')[0]
    n_replace = int(round(sb_frac * len(bi_idx)))
    np.random.seed(seed)
    chosen = np.random.choice(bi_idx, size=n_replace, replace=False)
    new_syms = list(syms)
    for idx in chosen:
        new_syms[idx] = 'Sb'
    supercell.set_chemical_symbols(new_syms)
    c = Counter(supercell.get_chemical_symbols())
    total_cation = c.get('Bi', 0) + c.get('Sb', 0)
    actual_pct = c.get('Sb', 0) / total_cation * 100 if total_cation > 0 else 0
    return supercell, actual_pct, c


def write_run_in(path, ec, transport_dir, nx, ny, nz,
                 L_transport, A_cross, n_atoms, mat_name, shc_group=4):
    t_axis   = DIR_MAP[transport_dir]
    pot_name = os.path.basename(POTENTIAL)
    txt = f"""potential {pot_name} 0

# Material    : {mat_name}
# Transport   : {transport_dir.upper()}
# Supercell   : {nx} x {ny} x {nz}
# L_transport : {L_transport/10:.3f} nm
# A_cross     : {A_cross/100:.4f} nm2
# Atoms       : {n_atoms}

minimize fire 1.0e-5 100000 1
ensemble    nve
time_step   0
dump_xyz    -1 0 1 relaxed.xyz
run         1

velocity {T_TARGET}
time_step 1
ensemble npt_ber {T_TARGET} {T_TARGET} {THERMO_COUPLING} 0.0 0.0 0.0 0.0 0.0 0.0 {ec['C11']:.1f} {ec['C22']:.1f} {ec['C33']:.1f} {ec['C44']:.1f} {ec['C55']:.1f} {ec['C66']:.1f} 1000
dump_thermo 1000
run {EQ_STEPS}

fix 0
ensemble heat_lan {T_TARGET} {THERMO_COUPLING} {T_DELTA} {G_SRC} {G_SNK}
compute 0 10 100 temperature jp jk
compute_shc 2 250 {t_axis} 1000 {MAX_OMEGA} group 0 {shc_group}
dump_thermo 1000
run {PROD_STEPS}
"""
    with open(path, 'w') as f:
        f.write(txt)


def write_setup_txt(path, mat_name, transport_dir, nx, ny, nz,
                    lx, ly, lz, L_transport, A_cross,
                    group_ids, z_bounds,
                    sb_frac=None, actual_sb_pct=None, comp=None):
    atoms_per_group = {g: int(np.sum(group_ids == g)) for g in range(9)}
    roles = {0:'wall (fixed)', 1:'src  (heat+)', 8:'snk  (heat-)'}
    for g in range(2, 8):
        roles[g] = f'mid  G{g}     '

    lines = ["=" * 65, f"NEMD SETUP  |  {mat_name}", "=" * 65]
    lines += [
        f"  Material       : {mat_name}",
        f"  Transport dir  : {transport_dir.upper()}",
    ]
    if sb_frac is not None:
        lines += [
            f"  Sb fraction    : {sb_frac*100:.0f}%  (actual {actual_sb_pct:.2f}%)",
            f"  Bi atoms       : {comp.get('Bi',0)}",
            f"  Sb atoms       : {comp.get('Sb',0)}",
            f"  Te atoms       : {comp.get('Te',0)}",
        ]
    lines += [
        f"  Supercell      : {nx} x {ny} x {nz}",
        f"  Total atoms    : {sum(atoms_per_group.values())}",
        f"  Lx             : {lx:.4f} A  ({lx/10:.4f} nm)",
        f"  Ly             : {ly:.4f} A  ({ly/10:.4f} nm)",
        f"  Lz             : {lz:.4f} A  ({lz/10:.4f} nm)",
        f"  L_transport    : {L_transport:.6f} A  ({L_transport/10:.6f} nm)",
        f"  A_cross        : {A_cross:.6f} A2  ({A_cross/100:.6f} nm2)",
        f"  T              : {T_TARGET} K  |  dT = {T_DELTA} K",
        f"  T_source       : {T_TARGET+T_DELTA} K  (group {G_SRC})",
        f"  T_sink         : {T_TARGET-T_DELTA} K  (group {G_SNK})",
        "",
        "GROUP BOUNDARIES  (along transport direction)",
        f"  {'GID':<5} {'Role':<16} {'z_min (A)':>12} {'z_max (A)':>12} "
        f"{'L (A)':>10} {'L (nm)':>8} {'Atoms':>7}",
        "  " + "-"*73,
    ]
    z0, z1 = z_bounds[0], z_bounds[1]
    lines.append(f"  {0:<5} {'wall-L (fixed)':<16} {z0:>12.4f} {z1:>12.4f} "
                 f"{z1-z0:>10.4f} {(z1-z0)/10:>8.4f} {atoms_per_group[0]//2:>7}")
    for seg in range(1, 9):
        z_lo, z_hi = z_bounds[seg], z_bounds[seg+1]
        lines.append(f"  {seg:<5} {roles[seg]:<16} {z_lo:>12.4f} {z_hi:>12.4f} "
                     f"{z_hi-z_lo:>10.4f} {(z_hi-z_lo)/10:>8.4f} {atoms_per_group[seg]:>7}")
    z_lo, z_hi = z_bounds[9], z_bounds[10]
    lines.append(f"  {0:<5} {'wall-R (fixed)':<16} {z_lo:>12.4f} {z_hi:>12.4f} "
                 f"{z_hi-z_lo:>10.4f} {(z_hi-z_lo)/10:>8.4f} {atoms_per_group[0]//2:>7}")

    z_src_cen = (z_bounds[1] + z_bounds[2]) / 2.0
    z_snk_cen = (z_bounds[9] + z_bounds[8]) / 2.0
    L_eff = z_snk_cen - z_src_cen
    lines += [
        "",
        "L_EFF  (src center to snk center, for kappa calculation)",
        f"  z_src_center   : {z_src_cen:.4f} A",
        f"  z_snk_center   : {z_snk_cen:.4f} A",
        f"  L_eff          : {L_eff:.4f} A  ({L_eff/10:.4f} nm)",
        "",
        "BIN CENTERS  (use for dT/dx gradient fit)",
        f"  {'GID':<5} {'Role':<16} {'center (A)':>12} {'center (nm)':>12} {'frac of L':>10}",
        "  " + "-"*57,
    ]
    for seg in range(1, 9):
        cen = (z_bounds[seg] + z_bounds[seg+1]) / 2.0
        lines.append(f"  {seg:<5} {roles[seg]:<16} {cen:>12.4f} "
                     f"{cen/10:>12.4f} {cen/L_transport:>10.4f}")
    lines.append("=" * 65)
    with open(path, 'w') as f:
        f.write("\n".join(lines))


def copy_files(dst_folder):
    for src in [POTENTIAL, SUBMIT_SH]:
        if src is None: continue
        fname = os.path.basename(src)
        dst   = os.path.join(dst_folder, fname)
        if os.path.exists(src) and not os.path.exists(dst):
            shutil.copy(src, dst)


# ============================================================
# MAIN — PURE MATERIALS
# ============================================================
print("=" * 70)
print("MASTER NEMD SETUP — TARGET A_CROSS = 120 nm² (ALL DIRECTIONS)")
print("=" * 70)

summary = []

pure_materials = [
    ('Bi2Te3', BI2TE3_PWI, EC_Bi2Te3),
    ('Sb2Te3', SB2TE3_PWI, EC_Sb2Te3),
]

for mat_name, pwi_file, ec in pure_materials:
    print(f"\n{'='*70}")
    print(f"BUILDING: {mat_name}")
    print(f"{'='*70}")

    atoms = read(pwi_file)
    uc    = atoms.cell.lengths()
    print(f"  Unit cell: a={uc[0]:.3f}  b={uc[1]:.3f}  c={uc[2]:.3f} A  "
          f"({len(atoms)} atoms)")

    for tdir in TRANSPORT_DIRS:
        taxis = DIR_MAP[tdir]
        print(f"\n  Direction: {tdir.upper()}")

        # Compute cross reps ONCE — locked to TARGET_A_CROSS = 120 nm²
        _, _, _, _, cross_reps = build_supercell(
            atoms, tdir, min(TARGET_LENGTHS), fixed_cross_reps=None)
        
        # Report what area we'll actually get
        _sc_test, _, _, _, _ = build_supercell(
            atoms, tdir, min(TARGET_LENGTHS), fixed_cross_reps=cross_reps)
        _A_test = get_A_cross(_sc_test, tdir) / 100.0
        print(f"    Cross reps locked: {cross_reps}  →  A_cross ≈ {_A_test:.3f} nm²  (target 120.0 nm²)")

        for target_L in TARGET_LENGTHS:
            folder_name = f"L_{int(target_L*10):04d}A"
            out_folder  = os.path.join(mat_name, tdir.upper(), folder_name)
            os.makedirs(out_folder, exist_ok=True)

            sc, nx, ny, nz, _ = build_supercell(
                atoms, tdir, target_L, fixed_cross_reps=cross_reps)

            lengths     = sc.cell.lengths()
            lx, ly, lz  = lengths
            L_transport = lengths[taxis]
            A_cross     = get_A_cross(sc, tdir)

            gids, zbounds, _ = assign_groups(sc, taxis)
            sc.arrays['group'] = gids

            write(os.path.join(out_folder, 'model.xyz'), sc, format='extxyz')
            write_run_in(os.path.join(out_folder, 'run.in'),
                         ec, tdir, nx, ny, nz, L_transport, A_cross,
                         len(sc), mat_name)
            write_setup_txt(os.path.join(out_folder, 'nemd_setup.txt'),
                            mat_name, tdir, nx, ny, nz,
                            lx, ly, lz, L_transport, A_cross, gids, zbounds)
            copy_files(out_folder)

            n_src = np.sum(gids == G_SRC)
            n_snk = np.sum(gids == G_SNK)
            ratio = n_src/n_snk if n_snk > 0 else 0
            flag  = "✅" if 0.8 < ratio < 1.2 else "⚠️"
            print(f"    {flag} {out_folder}  |  {len(sc):>7} atoms  |  "
                  f"L={L_transport/10:.3f} nm  |  A={A_cross/100:.3f} nm²")

            summary.append({'folder': out_folder, 'material': mat_name,
                            'dir': tdir.upper(), 'L_nm': L_transport/10,
                            'A_nm2': A_cross/100, 'n_atoms': len(sc), 'sb_pct': 0.0})


# ============================================================
# MAIN — ALLOYS
# ============================================================
print(f"\n{'='*70}")
print(f"BUILDING: ALLOYS")
print(f"{'='*70}")

for sb_frac in SB_FRACTIONS:
    sb_pct     = int(round(sb_frac * 100))
    alloy_name = f"BiSbTe{sb_pct}"
    ec         = interp_ec(sb_frac)
    print(f"\n  --- {alloy_name} ---")

    for tdir in TRANSPORT_DIRS:
        taxis = DIR_MAP[tdir]

        for target_L in TARGET_LENGTHS:
            folder_name = f"L_{int(target_L*10):04d}A"
            src_model   = os.path.join('Bi2Te3', tdir.upper(), folder_name, 'model.xyz')
            if not os.path.exists(src_model):
                print(f"    SKIP {alloy_name}/{tdir.upper()}/{folder_name} — source not found")
                continue

            out_folder = os.path.join(alloy_name, tdir.upper(), folder_name)
            os.makedirs(out_folder, exist_ok=True)

            sc          = read(src_model)
            lengths     = sc.cell.lengths()
            lx, ly, lz  = lengths
            L_transport = lengths[taxis]
            A_cross     = get_A_cross(sc, tdir)

            gids_orig = sc.arrays.get('group', None)
            if gids_orig is None:
                gids_orig, _, _ = assign_groups(sc, taxis)

            sc, actual_pct, comp = substitute_bi_sb(sc, sb_frac, RANDOM_SEED)
            sc.arrays['group'] = gids_orig
            zbounds = get_z_bounds(L_transport)

            uc_ref = read(BI2TE3_PWI).cell.lengths()
            nx = int(round(lx / uc_ref[0]))
            ny = int(round(ly / uc_ref[1]))
            nz = int(round(lz / uc_ref[2]))

            write(os.path.join(out_folder, 'model.xyz'), sc, format='extxyz')
            write_run_in(os.path.join(out_folder, 'run.in'),
                         ec, tdir, nx, ny, nz, L_transport, A_cross,
                         len(sc), alloy_name)
            write_setup_txt(os.path.join(out_folder, 'nemd_setup.txt'),
                            alloy_name, tdir, nx, ny, nz,
                            lx, ly, lz, L_transport, A_cross,
                            gids_orig, zbounds,
                            sb_frac=sb_frac, actual_sb_pct=actual_pct, comp=comp)
            copy_files(out_folder)

            with open(os.path.join(out_folder, 'stoichiometry.txt'), 'w') as f:
                f.write(f"Material    : {alloy_name}\nSb target   : {sb_pct}%\n"
                        f"Sb actual   : {actual_pct:.2f}%\nSeed        : {RANDOM_SEED}\n"
                        f"Bi atoms    : {comp.get('Bi',0)}\nSb atoms    : {comp.get('Sb',0)}\n"
                        f"Te atoms    : {comp.get('Te',0)}\nDirection   : {tdir.upper()}\n"
                        f"L_transport : {L_transport:.4f} A = {L_transport/10:.4f} nm\n"
                        f"A_cross     : {A_cross:.4f} A2 = {A_cross/100:.4f} nm2\n")

            print(f"    ✅ {out_folder}  |  {len(sc):>7} atoms  |  "
                  f"L={L_transport/10:.3f} nm  |  A={A_cross/100:.3f} nm²  |  Sb={actual_pct:.1f}%")

            summary.append({'folder': out_folder, 'material': alloy_name,
                            'dir': tdir.upper(), 'L_nm': L_transport/10,
                            'A_nm2': A_cross/100, 'n_atoms': len(sc), 'sb_pct': actual_pct})


# ============================================================
# CROSS-SECTION VERIFICATION
# ============================================================
print(f"\n{'='*70}")
print("CROSS-SECTION CONSISTENCY CHECK")
print(f"{'='*70}")
from itertools import groupby
key = lambda r: (r['material'], r['dir'])
for (mat, d), grp in groupby(sorted(summary, key=key), key=key):
    rows  = list(grp)
    areas = [r['A_nm2'] for r in rows]
    ok    = max(areas) - min(areas) < 0.01
    near120 = abs(areas[0] - 120.0) < 5.0
    flag  = "✅" if ok else "❌ INCONSISTENT!"
    flag2 = f"  {'~120 nm²' if near120 else '⚠️ FAR FROM 120'}"
    print(f"  {mat:12} {d}:  A = {areas[0]:.4f} nm²  (all {len(rows)} lengths)  {flag}{flag2}")

# ============================================================
# FINAL SUMMARY
# ============================================================
print(f"\n{'='*70}")
print(f"  {'Folder':<35} {'Dir':<4} {'L(nm)':<8} {'A(nm²)':<10} {'Atoms':<8} Sb%")
print("  " + "-"*70)
for r in summary:
    sb_str = f"{r['sb_pct']:.1f}%" if r['sb_pct'] > 0 else "pure"
    print(f"  {r['folder']:<35} {r['dir']:<4} {r['L_nm']:<8.3f} "
          f"{r['A_nm2']:<10.4f} {r['n_atoms']:<8} {sb_str}")
print(f"\n  Total folders: {len(summary)}")
print(f"{'='*70}")
