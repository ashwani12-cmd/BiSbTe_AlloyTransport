"""
add_grouping_method1.py
=======================
Adds a second grouping method (method 1) to a GPUMD model.xyz file
for use with compute_shc to get spectral interface conductance G(ω).

Your system layout (confirmed from file analysis):
───────────────────────────────────────────────────────────────
 g0   atoms    z-range (Å)     region
  0    1260    0   – 598       Wall atoms (fixed, both ends)
  1    3500    10  –  58       Heat source (Langevin hot, Bi₂Te₃ side)
  2    5600    60  – 139       Bi₂Te₃ bulk bin 1
  3    5600   141  – 219       Bi₂Te₃ bulk bin 2
  4    5600   222  – 300       Bi₂Te₃ INTERFACE BIN  ← SHC LEFT
  5    5600   302  – 379       Sb₂Te₃ INTERFACE BIN  ← SHC RIGHT
  6    5740   382  – 462       Sb₂Te₃ bulk bin 2
  7    5600   463  – 541       Sb₂Te₃ bulk bin 3
  8    3500   543  – 591       Heat sink (Langevin cold, Sb₂Te₃ side)
───────────────────────────────────────────────────────────────

Method 1 assignment:
  g0 = 4  →  g1 = 1   (Bi₂Te₃ interface, SHC LEFT)
  g0 = 5  →  g1 = 2   (Sb₂Te₃ interface, SHC RIGHT)
  everything else  →  g1 = 0   (excluded from SHC)

After running, update run.in:
  OLD: compute_shc 2 250 2 1000 50.0 group 0 4
  NEW: compute_shc 2 250 2 1000 50.0 group 1 -1

  group 1 -1  =  all method-1 groups EXCEPT g1=0
              =  g1=1 (Bi₂Te₃) + g1=2 (Sb₂Te₃)
              =  true cross-interface spectral conductance G(ω) ✅

Usage:
  python3 add_grouping_method1.py
  (run in same folder as model.xyz)
"""

import os
import sys
import re
import shutil
from collections import defaultdict

# ══════════════════════════════════════════════════════════════
# USER SETTINGS  — change these if your group layout is different
# ══════════════════════════════════════════════════════════════
INPUT_FILE  = 'model.xyz'
OUTPUT_FILE = 'model.xyz'          # overwrite in-place
BACKUP_FILE = 'model_backup.xyz'   # original saved here first

# Method 0 group ID  →  Method 1 group ID
GROUP_MAP = {
    0: 0,   # wall              → excluded
    1: 0,   # heat source       → excluded
    2: 0,   # Bi₂Te₃ bulk bin1 → excluded
    3: 0,   # Bi₂Te₃ bulk bin2 → excluded
    4: 1,   # Bi₂Te₃ interface → SHC LEFT  ✅
    5: 2,   # Sb₂Te₃ interface → SHC RIGHT ✅
    6: 0,   # Sb₂Te₃ bulk bin2 → excluded
    7: 0,   # Sb₂Te₃ bulk bin3 → excluded
    8: 0,   # heat sink         → excluded
}
# ══════════════════════════════════════════════════════════════


def parse_model(path):
    """Read model.xyz, return n_atoms, header, list of atom dicts."""
    with open(path, 'r') as f:
        lines = f.readlines()

    n_atoms = int(lines[0].strip())
    header  = lines[1].rstrip('\n')

    atoms = []
    for i, line in enumerate(lines[2:2 + n_atoms]):
        p = line.split()
        if len(p) < 5:
            print(f"  WARNING: skipping malformed line {i+3}: {line.rstrip()}")
            continue
        atoms.append({
            'species': p[0],
            'x':       float(p[1]),
            'y':       float(p[2]),
            'z':       float(p[3]),
            'g0':      int(p[4]),
        })

    return n_atoms, header, atoms


def update_header(header):
    """Change group:I:1 → group:I:2 (or verify already :2)."""
    if 'group:I:1' in header:
        return header.replace('group:I:1', 'group:I:2'), 'upgraded group:I:1 → group:I:2'
    elif 'group:I:2' in header:
        return header, 'already group:I:2 — will overwrite method-1 column'
    else:
        return None, 'ERROR: no group:I column found in header'


def assign_method1(g0, group_map):
    """Return method-1 group id for a given method-0 group id."""
    if g0 not in group_map:
        print(f"  WARNING: unknown g0={g0} — assigning g1=0 (excluded)")
        return 0
    return group_map[g0]


def write_model(path, n_atoms, header, atoms):
    """Write updated model.xyz with two group columns."""
    with open(path, 'w') as f:
        f.write(f"{n_atoms}\n")
        f.write(f"{header}\n")
        for a in atoms:
            f.write(
                f"{a['species']:<4s}  "
                f"{a['x']:20.14f}  "
                f"{a['y']:20.14f}  "
                f"{a['z']:20.14f}  "
                f"{a['g0']:4d}  {a['g1']:4d}\n"
            )


def print_summary(atoms, group_map):
    """Print per-group summary table."""
    g0_info = defaultdict(lambda: {'count': 0, 'g1': 0,
                                   'zmin': 1e9, 'zmax': -1e9,
                                   'species': defaultdict(int)})
    g1_count = defaultdict(int)

    for a in atoms:
        g = a['g0']
        g0_info[g]['count']  += 1
        g0_info[g]['g1']      = a['g1']
        g0_info[g]['species'][a['species']] += 1
        if a['z'] < g0_info[g]['zmin']: g0_info[g]['zmin'] = a['z']
        if a['z'] > g0_info[g]['zmax']: g0_info[g]['zmax'] = a['z']
        g1_count[a['g1']] += 1

    region_names = {
        0: 'Wall (fixed)',
        1: 'Heat source',
        2: 'Bi₂Te₃ bulk bin 1',
        3: 'Bi₂Te₃ bulk bin 2',
        4: 'Bi₂Te₃ interface  ← SHC LEFT',
        5: 'Sb₂Te₃ interface  ← SHC RIGHT',
        6: 'Sb₂Te₃ bulk bin 2',
        7: 'Sb₂Te₃ bulk bin 3',
        8: 'Heat sink',
    }

    W = 75
    print(f"\n{'='*W}")
    print(f"  Method 0 groups  →  Method 1 assignment")
    print(f"{'='*W}")
    print(f"  {'g0':>3}  {'atoms':>7}  {'g1':>4}  {'z_min':>8}  {'z_max':>8}  region")
    print(f"  {'-'*65}")
    for g in sorted(g0_info):
        v    = g0_info[g]
        name = region_names.get(g, f'unknown group {g}')
        flag = ' ✅' if v['g1'] != 0 else ''
        print(f"  {g:3d}  {v['count']:7d}  {v['g1']:4d}  "
              f"{v['zmin']:8.2f}  {v['zmax']:8.2f}  {name}{flag}")

    print(f"\n{'='*W}")
    print(f"  Method 1 summary  (used by compute_shc group 1 -1)")
    print(f"{'='*W}")
    print(f"  g1=0  {g1_count[0]:6d} atoms  → excluded from SHC")
    print(f"  g1=1  {g1_count[1]:6d} atoms  → Bi₂Te₃ interface (SHC LEFT)")
    print(f"  g1=2  {g1_count[2]:6d} atoms  → Sb₂Te₃ interface (SHC RIGHT)")
    print(f"  {'─'*40}")
    print(f"  Total {sum(g1_count.values()):6d} atoms")
    print(f"{'='*W}")


def main():
    # ── check input exists ───────────────────────────────────
    if not os.path.exists(INPUT_FILE):
        print(f"ERROR: '{INPUT_FILE}' not found!")
        print("Run this script in the same directory as your model.xyz")
        sys.exit(1)

    # ── backup original ──────────────────────────────────────
    shutil.copy2(INPUT_FILE, BACKUP_FILE)
    print(f"Backup saved → {BACKUP_FILE}")

    # ── parse ────────────────────────────────────────────────
    print(f"Reading {INPUT_FILE}...")
    n_atoms, header, atoms = parse_model(INPUT_FILE)
    print(f"  Atoms   : {n_atoms}")
    print(f"  Groups  : {sorted(set(a['g0'] for a in atoms))}")

    # ── update header ────────────────────────────────────────
    new_header, msg = update_header(header)
    print(f"  Header  : {msg}")
    if new_header is None:
        print("  ABORTING — cannot proceed without group column")
        sys.exit(1)

    # ── assign method 1 ──────────────────────────────────────
    for a in atoms:
        a['g1'] = assign_method1(a['g0'], GROUP_MAP)

    # ── write output ─────────────────────────────────────────
    write_model(OUTPUT_FILE, n_atoms, new_header, atoms)
    print(f"  Written → {OUTPUT_FILE}")

    # ── verify atom count ────────────────────────────────────
    with open(OUTPUT_FILE) as f:
        written = sum(1 for _ in f) - 2   # subtract header lines
    if written == n_atoms:
        print(f"  Atom count: {written} ✅")
    else:
        print(f"  Atom count MISMATCH: wrote {written}, expected {n_atoms} ❌")
        sys.exit(1)

    # ── summary table ────────────────────────────────────────
    print_summary(atoms, GROUP_MAP)

    # ── final instructions ───────────────────────────────────
    print(f"\n  ✅ Done! Now update run.in:")
    print(f"")
    print(f"     OLD:  compute_shc 2 250 2 1000 50.0 group 0 4")
    print(f"     NEW:  compute_shc 2 250 2 1000 50.0 group 1 -1")
    print(f"")
    print(f"  group 1 -1  =  method-1 groups 1 + 2")
    print(f"              =  Bi₂Te₃ interface + Sb₂Te₃ interface")
    print(f"              =  true spectral interface conductance G(ω) ✅")


if __name__ == '__main__':
    main()
