"""
Single-Interface Bi2Te3/Sb2Te3 NEMD Structure Builder
======================================================
For Kapitza resistance measurement you need EXACTLY ONE interface
in the measurement region. This script builds:

  |wall|src(Bi2Te3)|pure Bi2Te3|INTERFACE|pure Sb2Te3|snk(Sb2Te3)|wall|

Structure along z:
  - Left half  : pure Bi2Te3  (~50% of cell)
  - Right half : pure Sb2Te3  (~50% of cell)
  - 1 interface at z ~ 0.50 c
  - Periodic image interface at z=0/c is SUPPRESSED by wall groups

Group assignment:
  G0 (wall) : thin slabs at both ends — fixed, no thermostat
  G1 (src)  : inside Bi2Te3 region   — heated to T+dT
  G2-G7     : mid bins across both materials + interface
  G8 (snk)  : inside Sb2Te3 region   — cooled to T-dT

The periodic interface (at z=0) is inside G0 (wall/fixed) so it
does NOT contribute to the measured ΔT — only the real interface
at z~0.5c is measured.

Usage:
    python3 build_single_interface.py
    # Generates structures + run.in + nemd_setup.txt for each config
"""

import os
import shutil
import glob
import numpy as np
from ase import Atoms
from ase.io import read, write

# ============================================================
# SETTINGS
# ============================================================
BI2TE3_XYZ = "./NEMD_supercells/4QL_TeTe_10x10x10.xyz"   # use as template
SB2TE3_XYZ = None   # will be built from Sb half of TeTe structure

# Actually we build from the original conv cells
BI2TE3_PWI = "./espresso_Bi2Te3_conventional.pwi"
SB2TE3_PWI = "./espresso_Sb2Te3_conventional.pwi"

OUT_DIR = "./superlattice_NEMD_singleIF_FIXED"

# Auto-detect NEP
_nep = glob.glob("./nep*.txt")
if not _nep:
    raise FileNotFoundError("No nep*.txt found!")
POTENTIAL = sorted(_nep)[-1]
print(f"  NEP : {POTENTIAL}")

_sub = glob.glob("./submit.sh") + glob.glob("./*.sh")
SUBMIT_SH = _sub[0] if _sub else None

# MD settings
T_TARGET        = 300
T_DELTA         = 30    # src = T+30, snk = T-30  => 60 K total across system (Chowdhury 2021)
THERMO_COUPLING = 100
MAX_OMEGA       = 50.0
EQ_STEPS        = 1_000_000
PROD_STEPS      = 5_000_000
EC = dict(C11=55.0, C22=55.0, C33=10.0, C44=8.0, C55=8.0, C66=20.0)

# Group fractions
# Wall at both ends suppresses periodic interface
# Src deep in Bi2Te3, Snk deep in Sb2Te3
FRAC_WALL = 0.05   # 5% each end — buries periodic interface at z=0
FRAC_SRC  = 0.15   # src thermostat in Bi2Te3 (z ~ 0.05-0.20)
FRAC_SNK  = 0.15   # snk thermostat in Sb2Te3 (z ~ 0.80-0.95)
# Group scheme: FIXED QL counts at each end (paper: 1 QL fixed + 5 QL bath),
# so across the length series only the middle measurement region grows.
N_FIX_QL  = 1      # fixed wall (each end)  -- buries periodic interface at z=0
N_BATH_QL = 5      # Langevin bath (each end) -- >=5 QL for convergence (Chowdhury 2021)
N_MID     = 6      # measurement bins between the two baths (interface sits inside)
G_SRC = 1
G_SNK = 8

# Target sizes: vary n_bi and n_sb (number of QLs each side)
# In-plane: fixed 10x10 repeats
# For finite-size correction: vary total length (n_bi + n_sb)
# Terminations to compare: Te-Te and Metal-Metal
# Change this line in CONFIGS:
CONFIGS = [
    (10, 10, "TeTe", 10, 14),   # Lz ~20 nm
    (15, 15, "TeTe", 10, 14),   # Lz ~30 nm
    (30, 30, "TeTe", 10, 14),   # Lz ~60 nm  (longest, for 1/L extrapolation)
]

# ============================================================
# STRUCTURE BUILDERS
# ============================================================

def scale_to_a(conv, a_target):
    atoms = conv.copy()
    s = a_target / conv.cell[0, 0]
    c = atoms.cell.copy()
    c[0,0] *= s; c[1,0] *= s; c[1,1] *= s
    atoms.set_cell(c, scale_atoms=True)
    return atoms

def vdw_offset(conv_scaled):
    """z of the plane that STARTS a quintuple layer (the Te just after the
    widest interlayer gap = the van der Waals gap). Anchoring QL extraction
    here guarantees a symmetric Te-capped QL (Te-X-Te-X-Te)."""
    c   = conv_scaled.cell[2,2]
    pos = conv_scaled.positions.copy(); pos[:,2] %= c
    z   = np.sort(np.unique(np.round(pos[:,2], 3)))
    zc  = np.append(z, z[0] + c)
    gaps = np.diff(zc)
    i   = int(np.argmax(gaps))          # widest gap is AFTER plane z[i]
    return z[(i + 1) % len(z)]          # QL begins at the plane after the gap

def extract_QL(conv_scaled, ql_index=0):
    QL_t = conv_scaled.cell[2,2] / 3
    z0   = vdw_offset(conv_scaled)       # <-- FIX: anchor on real vdW gap, not z=0
    pos  = conv_scaled.positions.copy()
    syms = np.array(conv_scaled.get_chemical_symbols())
    c    = conv_scaled.cell[2,2]
    zz   = (pos[:,2] - z0) % c            # shift frame so QL start -> 0
    zz   = np.where(zz > c - 1e-3, zz - c, zz)   # snap anchor plane (~c from fp) back to ~0
    pos[:,2] = zz
    z_lo = ql_index * QL_t
    mask = (pos[:,2] >= z_lo - 0.05) & (pos[:,2] < z_lo + QL_t - 0.05)
    p = pos[mask].copy(); s = syms[mask].tolist()
    p[:,2] -= z_lo
    order = np.argsort(p[:,2])
    return [s[i] for i in order], p[order], QL_t

def build_slab(conv, n_QL, a_target, reverse=False):
    """Build a slab of n_QL quintuple layers."""
    scaled = scale_to_a(conv, a_target)
    all_syms=[]; all_pos=[]; z_off=0.0
    for _ in range(n_QL):
        syms, pos, QL_t = extract_QL(scaled, 0)
        p = pos.copy()
        if reverse:
            p[:,2] = QL_t - p[:,2]
            p[:,2] -= p[:,2].min()
            order = np.argsort(p[:,2])
            p = p[order]; syms = [syms[i] for i in order]
        p[:,2] += z_off
        all_pos.append(p); all_syms += syms; z_off += QL_t
    return all_syms, np.vstack(all_pos), z_off

def build_single_interface(n_bi, n_sb, termination, na, nb, bi_conv, sb_conv):
    """
    Build a single-interface Bi2Te3/Sb2Te3 slab.
    termination='TeTe'       → Te | Te at interface (natural)
    termination='MetalMetal' → Bi | Sb at interface
    In-plane: na x nb repeats.
    """
    a_avg = (bi_conv.cell[0,0] + sb_conv.cell[0,0]) / 2

    # Build Bi2Te3 slab (normal orientation — Te on top)
    bi_syms, bi_pos, bi_c = build_slab(bi_conv, n_bi, a_avg, reverse=False)

    # Build Sb2Te3 slab
    # For TeTe: normal → bottom layer is Te, top layer is Te
    #           Bi2Te3 top (Te) meets Sb2Te3 bottom (Te) → Te-Te ✅
    # For MetalMetal: reverse Sb2Te3 → bottom layer becomes Sb
    #           Bi2Te3 top (Bi) meets Sb2Te3 bottom (Sb) → Metal-Metal ✅
    reverse_sb = (termination == "MetalMetal")
    sb_syms, sb_pos, sb_c = build_slab(sb_conv, n_sb, a_avg, reverse=reverse_sb)
    sb_pos[:,2] += bi_c

    all_syms = bi_syms + sb_syms
    all_pos  = np.vstack([bi_pos, sb_pos])
    c_total  = bi_c + sb_c

    # Build 1x1 cell
    s = a_avg / bi_conv.cell[0,0]
    cell = bi_conv.cell.copy()
    cell[0,0] *= s; cell[1,0] *= s; cell[1,1] *= s
    cell[2,2]  = c_total

    unit = Atoms(symbols=all_syms, positions=all_pos, cell=cell, pbc=True)

    # Repeat in-plane
    big = unit.repeat([na, nb, 1])
    return big, bi_c, sb_c, c_total

def assign_groups_single_IF(atoms, bi_c, c_total, n_bi, n_sb):
    """
    Fixed-QL groups for single-interface NEMD (paper-style):
      [1 QL wall][5 QL hot bath][ ...N_MID measurement bins... ][5 QL cold bath][1 QL wall]
        G0            G1                 G2..G7                       G8            G0
    Bottom wall+bath are in Bi2Te3 (QL_t_Bi), top wall+bath in Sb2Te3 (QL_t_Sb).
    Only the middle measurement region grows with system length; reservoirs stay fixed.
    Periodic interface at z=0 is buried inside the 1-QL walls.
    """
    z    = atoms.positions[:,2]
    c    = atoms.cell[2,2]
    QL_Bi = bi_c / n_bi
    QL_Sb = (c_total - bi_c) / n_sb
    # nudge reservoir boundaries into the vdW gap (~0.13*QL below a QL edge)
    # so walls/baths contain whole QLs, never split an atom plane
    dBi = 0.13 * QL_Bi
    dSb = 0.13 * QL_Sb

    z_wallL = N_FIX_QL * QL_Bi - dBi                       # end of bottom wall
    z_srcR  = (N_FIX_QL + N_BATH_QL) * QL_Bi - dBi         # end of hot bath (G1)
    z_wallR = c - N_FIX_QL * QL_Sb + dSb                   # start of top wall
    z_snkL  = c - (N_FIX_QL + N_BATH_QL) * QL_Sb + dSb     # start of cold bath (G8)

    # N_MID measurement bins between the two baths
    mid_edges = np.linspace(z_srcR, z_snkL, N_MID + 1)
    zb = [0.0, z_wallL, z_srcR] + list(mid_edges[1:-1]) + [z_snkL, z_wallR, c]
    # zb has 11 entries: [0, wallL, srcR, m1..m5, snkL, wallR, c]

    gids = np.zeros(len(atoms), dtype=int)
    for i, zi in enumerate(z):
        if   zi < zb[1] or zi >= zb[9]:   gids[i] = 0     # walls
        elif zi < zb[2]:                   gids[i] = 1     # hot bath
        elif zi < zb[3]:                   gids[i] = 2
        elif zi < zb[4]:                   gids[i] = 3
        elif zi < zb[5]:                   gids[i] = 4
        elif zi < zb[6]:                   gids[i] = 5
        elif zi < zb[7]:                   gids[i] = 6
        elif zi < zb[8]:                   gids[i] = 7
        else:                              gids[i] = 8     # cold bath

    if_frac = bi_c / c_total
    if_z    = bi_c
    if_group = None
    for g in range(2, 8):
        if zb[g] <= if_z < zb[g+1]:
            if_group = g
            break
    bf = [b / c for b in zb]
    return gids, zb, bf, if_frac, if_group

def write_run_in(path, label, nx, ny, nz_dummy,
                 L_transport, A_cross, n_atoms):
    pot_name = os.path.basename(POTENTIAL)
    txt = f"""potential {pot_name} 0

# Superlattice  : {label}
# L_transport   : {L_transport/10:.3f} nm  (z-direction)
# A_cross       : {A_cross/100:.4f} nm²
# Atoms         : {n_atoms}
# Interfaces    : 1 (real) + 1 (periodic, buried in wall)

# ── Step 1: Minimize ─────────────────────────────────────────
minimize fire 1.0e-5 100000 1
ensemble    nve
time_step   0
dump_xyz    -1 0 1 relaxed.xyz
run         1

# ── Step 2: NPT Equilibration ─────────────────────────────────
velocity {T_TARGET}
time_step 1
ensemble npt_ber {T_TARGET} {T_TARGET} {THERMO_COUPLING} 0.0 0.0 0.0 0.0 0.0 0.0 {EC['C11']:.1f} {EC['C22']:.1f} {EC['C33']:.1f} {EC['C44']:.1f} {EC['C55']:.1f} {EC['C66']:.1f} 1000
dump_thermo 1000
run {EQ_STEPS}

# ── Step 3: NEMD Production ───────────────────────────────────
fix 0
ensemble heat_lan {T_TARGET} {THERMO_COUPLING} {T_DELTA} {G_SRC} {G_SNK}
compute 0 10 100 temperature jp jk
compute_shc 2 250 2 1000 {MAX_OMEGA} group 0 4
dump_thermo 1000
run {PROD_STEPS}
"""
    with open(path, 'w') as f:
        f.write(txt)

def write_setup_txt(path, label, n_bi, n_sb, na, nb,
                    lx, ly, lz, A_cross, gids, zb, bf,
                    if_frac, if_group, bi_c, sb_c):
    apg = {g: int(np.sum(gids == g)) for g in range(9)}
    roles = {0:'wall (fixed)', 1:'src (Bi2Te3)', 8:'snk (Sb2Te3)'}
    for g in range(2,8): roles[g] = f'mid G{g}'

    z_src_cen = (zb[1] + zb[2]) / 2
    z_snk_cen = (zb[8] + zb[9]) / 2
    L_eff = z_snk_cen - z_src_cen

    lines = [
        "="*65,
        f"SINGLE-INTERFACE NEMD  |  {label}",
        "="*65,
        f"  Bi2Te3 slab   : {n_bi} QLs  ({bi_c:.2f} Å = {bi_c/10:.2f} nm)",
        f"  Sb2Te3 slab   : {n_sb} QLs  ({sb_c:.2f} Å = {sb_c/10:.2f} nm)",
        f"  In-plane reps : {na} x {nb}",
        f"  Total atoms   : {len(gids):,}",
        f"  Lz (total)    : {lz:.2f} Å  ({lz/10:.2f} nm)",
        f"  A_cross       : {A_cross:.2f} Å²  ({A_cross/100:.3f} nm²)",
        f"  T             : {T_TARGET} K  |  dT = {T_DELTA} K",
        f"  Src (G1)      : {T_TARGET+T_DELTA} K  in Bi2Te3",
        f"  Snk (G8)      : {T_TARGET-T_DELTA} K  in Sb2Te3",
        "",
        f"  INTERFACE     : z = {if_frac*lz:.2f} Å  (z/c = {if_frac:.4f})",
        f"  Interface in  : Group {if_group}  ← ΔT discontinuity here",
        f"  Periodic IF   : z = 0 Å  → buried in G0 wall ✅",
        "",
        "GROUP BOUNDARIES",
        f"  {'GID':<5} {'Role':<16} {'z_min':>10} {'z_max':>10} {'L(Å)':>8} {'L(nm)':>7} {'Atoms':>7}",
        "  "+"-"*65,
    ]
    lines.append(f"  {0:<5} {'wall-L':<16} {zb[0]:>10.2f} {zb[1]:>10.2f} "
                 f"{zb[1]-zb[0]:>8.2f} {(zb[1]-zb[0])/10:>7.3f} {apg[0]//2:>7}")
    for seg in range(1,9):
        z_lo, z_hi = zb[seg], zb[seg+1]
        iface = f"← IF @ {if_frac*lz:.1f}Å" if seg == if_group else ""
        lines.append(f"  {seg:<5} {roles[seg]:<16} {z_lo:>10.2f} {z_hi:>10.2f} "
                     f"{z_hi-z_lo:>8.2f} {(z_hi-z_lo)/10:>7.3f} {apg[seg]:>7}  {iface}")
    lines.append(f"  {0:<5} {'wall-R':<16} {zb[9]:>10.2f} {zb[10]:>10.2f} "
                 f"{zb[10]-zb[9]:>8.2f} {(zb[10]-zb[9])/10:>7.3f} {apg[0]//2:>7}")
    lines += [
        "",
        f"  L_eff (src→snk centers): {L_eff:.2f} Å  ({L_eff/10:.3f} nm)",
        "="*65,
    ]
    with open(path, 'w') as f:
        f.write("\n".join(lines))

# ============================================================
# MAIN
# ============================================================
bi_conv = read(BI2TE3_PWI)
sb_conv = read(SB2TE3_PWI)
a_avg   = (bi_conv.cell[0,0] + sb_conv.cell[0,0]) / 2

os.makedirs(OUT_DIR, exist_ok=True)
print(f"\n{'='*70}")
print("BUILDING SINGLE-INTERFACE STRUCTURES")
print(f"{'='*70}\n")

summary = []

for (n_bi, n_sb, term, na, nb) in CONFIGS:
    label = f"SI_{n_bi}Bi_{n_sb}Sb_{term}_{na}x{nb}"
    print(f"  Building: {label}")

    atoms, bi_c, sb_c, c_total = build_single_interface(
        n_bi, n_sb, term, na, nb, bi_conv, sb_conv)

    lx, ly, lz = atoms.cell.lengths()
    cell        = atoms.cell[:]
    A_cross     = np.linalg.norm(np.cross(cell[0], cell[1]))

    gids, zb, bf, if_frac, if_group = assign_groups_single_IF(
        atoms, bi_c, c_total, n_bi, n_sb)
    atoms.arrays['group'] = gids

    n_src = np.sum(gids == G_SRC)
    n_snk = np.sum(gids == G_SNK)
    ratio = n_src/n_snk if n_snk > 0 else 0

    # Sanity: interface must be in mid group (2-7), not in wall/src/snk
    if_in_mid = if_group is not None and 2 <= if_group <= 7
    # Periodic interface at z=0 must be in wall (G0)
    periodic_if_buried = bf[0] <= 0.0 < bf[1]  # always true by construction

    flag1 = "✅" if if_in_mid else "⚠️  IF NOT IN MID!"
    flag2 = "✅" if 0.8 < ratio < 1.2 else "⚠️  IMBALANCED"

    out = os.path.join(OUT_DIR, label)
    os.makedirs(out, exist_ok=True)

    write(os.path.join(out, "model.xyz"), atoms, format="extxyz")
    write_run_in(os.path.join(out, "run.in"),
                 label, na, nb, 1, lz, A_cross, len(atoms))
    write_setup_txt(os.path.join(out, "nemd_setup.txt"),
                    label, n_bi, n_sb, na, nb,
                    lx, ly, lz, A_cross, gids, zb, bf,
                    if_frac, if_group, bi_c, sb_c)

    for src in [POTENTIAL, SUBMIT_SH]:
        if src and os.path.exists(src):
            dst = os.path.join(out, os.path.basename(src))
            if not os.path.exists(dst):
                shutil.copy(src, dst)

    print(f"    {flag1} {flag2}")
    print(f"    Atoms     : {len(atoms):,}")
    print(f"    Lz        : {lz:.1f} Å ({lz/10:.1f} nm)  "
          f"[Bi={bi_c/10:.1f}nm | Sb={sb_c/10:.1f}nm]")
    print(f"    A_cross   : {A_cross/100:.3f} nm²")
    print(f"    Interface : z={if_frac*lz:.1f} Å (z/c={if_frac:.3f}) in G{if_group}")
    print(f"    Src/Snk   : {n_src}/{n_snk} atoms  ratio={ratio:.2f}")
    print()

    summary.append(dict(label=label, atoms=len(atoms),
                        Lz=lz/10, A=A_cross/100,
                        if_group=if_group, src=n_src, snk=n_snk))

print("="*70)
print(f"  {'Label':<40} {'Atoms':>8} {'Lz(nm)':>8} {'A(nm²)':>8} {'IF_grp':>7}")
print("  "+"-"*70)
for r in summary:
    print(f"  {r['label']:<40} {r['atoms']:>8,} {r['Lz']:>8.1f} "
          f"{r['A']:>8.3f} {r['if_group']:>7}")
print(f"\n  Total : {len(summary)} configs")
print(f"  Output: {OUT_DIR}/")
print("="*70)
