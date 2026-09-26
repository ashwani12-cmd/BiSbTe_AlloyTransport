#!/usr/bin/env bash
set -euo pipefail

rm -rf input_struct lmp_force_output lmp_format qe_forces
# ---------------- USER CONFIG ----------------
export LMP_EXEC="/home/IITB/multiscale-mechanics/ashwani12/software/builds/NEP_CPU/lammps-24Mar2022/src/lmp_mpi"

DIM="3 3 2"                      # phonopy supercell dim
QE_CELL="espresso_frac.pwi"       # QE cell in fractional coords
HEADER_FILE="header.in"           # LAMMPS header (pair_style, etc.)
INPUT_DIR="input_struct"
OUTPUT_DIR="lmp_format"
POT_DIR="pot_files"
FORCE_STAGE_DIR="lmp_force_output"
QE_FORCES_DIR="qe_forces"

# phonopy for FC2 + band
CONDA_PATH="/home/IITB/multiscale-mechanics/ashwani12/anaconda3/etc/profile.d/conda.sh"
CONDA_ENV="phonopy"

BAND_POINTS=51
# ---------------------------------------------

echo ">>> Cleaning old folders…"
rm -rf "${INPUT_DIR}" "${OUTPUT_DIR}" "${FORCE_STAGE_DIR}" "${QE_FORCES_DIR}"
mkdir -p "${INPUT_DIR}" "${OUTPUT_DIR}" "${POT_DIR}" "${FORCE_STAGE_DIR}" "${QE_FORCES_DIR}"

# sanity checks
[[ -f "${QE_CELL}" ]]     || { echo "ERROR: ${QE_CELL} not found in $(pwd)"; exit 1; }
[[ -f "${HEADER_FILE}" ]] || { echo "ERROR: ${HEADER_FILE} not found in $(pwd)"; exit 1; }

cp -n "${QE_CELL}" "${INPUT_DIR}/" || true
cp -n "${HEADER_FILE}" "${INPUT_DIR}/" || true

echo "====================================================="
echo "  STEP 1: phonopy displacements (QE, DIM=${DIM})"
echo "====================================================="

pushd "${INPUT_DIR}" >/dev/null

phonopy --qe -d --dim="${DIM}" -c "${QE_CELL}"

echo ">>> Generated phonopy_disp.yaml and supercell-*.in:"
ls supercell-*.in || { echo "ERROR: no supercell-*.in from phonopy"; popd >/dev/null; exit 1; }

mapfile -t IDX_ARR < <(ls supercell-*.in 2>/dev/null \
  | grep -Eo 'supercell-[0-9]{3,5}\.in' \
  | sed -E 's/.*-([0-9]{3,5})\.in/\1/' \
  | sort -n)

if (( ${#IDX_ARR[@]} == 0 )); then
  echo "ERROR: No supercell-XXXXX.in files found after phonopy."
  popd >/dev/null; exit 1
fi

MIN_IDX="${IDX_ARR[0]}"; MAX_IDX="${IDX_ARR[-1]}"
printf ">>> Found supercell range: %s … %s\n" "${MIN_IDX}" "${MAX_IDX}"

echo ">>> Building BiSbTe-XXXXX.in (QE input) by header+supercell…"
for n in $(seq $((10#${MIN_IDX})) $((10#${MAX_IDX}))); do
  pad=$(printf "%03d" "${n}")
  sc="supercell-${pad}.in"
  out="BiSbTe-${pad}.in"           # ← changed from Mg3Bi2
  if [[ -f "${sc}" ]]; then
    cat "${HEADER_FILE}" "${sc}" >| "${out}"
    echo "  ✓ ${out}"
  else
    echo "  ⚠️ Missing ${sc} — skipped"
  fi
done

popd >/dev/null

echo "====================================================="
echo "  STEP 2: QE → LAMMPS .data (ASE, Bi=1 Sb=2 Te=3)"
echo "====================================================="

echo ">>> Clearing old .data in ${OUTPUT_DIR}…"
find "${OUTPUT_DIR}" -maxdepth 1 -type f -name "*.data" -print -delete | sed 's/^/  🗑 /' || true

python3 - <<'PY'
import os, re, sys
from ase.io import read, write

INPUT_DIR  = "input_struct"
OUTPUT_DIR = "lmp_format"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── change species here for your composition ──
ALLOWED   = {"Bi", "Sb", "Te"}
SPECORDER = ["Bi", "Sb", "Te"]   # Bi=1, Sb=2, Te=3

pat = re.compile(r'^BiSbTe-(\d+)\.in$')   # ← changed pattern

files   = []
indices = []
for f in os.listdir(INPUT_DIR):
    m = pat.match(f)
    if m:
        files.append(f)
        indices.append(int(m.group(1)))

files = sorted(files)
if not files:
    sys.exit("No BiSbTe-XXX.in in input_struct.")

start = min(indices)
end   = max(indices)
print(f"  Range detected: BiSbTe-{start:03d}.in … BiSbTe-{end:03d}.in")

conv, skip = 0, 0
for idx in range(start, end+1):
    pad   = f"{idx:03d}"
    fname = f"BiSbTe-{pad}.in"
    src   = os.path.join(INPUT_DIR, fname)
    if not os.path.isfile(src):
        print(f"  ⚠️ Missing {fname} — skip")
        skip += 1
        continue

    dst = os.path.join(OUTPUT_DIR, f"BiSbTe-{pad}.data")
    try:
        atoms   = read(src, format="espresso-in")
        syms    = set(atoms.get_chemical_symbols())
        unknown = syms - ALLOWED
        if unknown:
            raise RuntimeError(
                f"Unsupported species in {fname}: {sorted(unknown)} "
                f"(allowed: {sorted(ALLOWED)})"
            )
        write(dst, atoms, format="lammps-data",
              atom_style="atomic", specorder=SPECORDER)
        print(f"  ✅ {fname} → {os.path.basename(dst)}  (Bi=1, Sb=2, Te=3)")
        conv += 1
    except Exception as e:
        print(f"  ❌ Error {fname}: {e}")

print(f"\nSummary: converted={conv}, skipped={skip}, out_dir={OUTPUT_DIR}")
PY


echo "====================================================="
echo "  STEP 3: create LAMMPS input files per .data"
echo "====================================================="

pushd "${OUTPUT_DIR}" >/dev/null
python3 - <<'PY'
import os, fnmatch, re

# ── NEP potential file ──
NEP_FILE = "nep_y2026_m04_d08_h20_m24_s14_generation100000.txt"   # ← updated

TEMPLATE = """units metal

box           tilt large
read_data supercell

mass 1 208.9804
mass 2 121.760
mass 3 127.600

pair_style nep
pair_coeff  * * ../pot_files/{nep}  Bi Sb Te

dump phonopy all custom 1 force.ID id type x y z fx fy fz
dump_modify phonopy format line "%d %d %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f" sort id
run 0
""".format(nep=NEP_FILE)

for f in list(os.listdir(".")):
    if os.path.isfile(f) and not fnmatch.fnmatch(f, "BiSbTe-*.data") and not f.endswith(".in"):
        try: os.remove(f)
        except: pass

data_files = sorted(fnmatch.filter(os.listdir("."), "BiSbTe-*.data"))
for data_file in data_files:
    m = re.search(r'BiSbTe-(\d+)\.data', data_file)
    if not m: continue
    idx      = m.group(1)
    in_name  = f"BiSbTe_in_{idx}.in"
    dump_name = f"force.{idx}"
    lines = []
    for line in TEMPLATE.splitlines(keepends=True):
        if line.strip().startswith("read_data"):
            lines.append(f"read_data {data_file}\n")
        elif line.strip().startswith("dump phonopy"):
            lines.append(f'dump phonopy all custom 1 {dump_name} id type x y z fx fy fz\n')
        else:
            lines.append(line.replace("force.ID", dump_name))
    with open(in_name, "w") as f:
        f.writelines(lines)
    print(f"  📝 {in_name} for {data_file}")
print("✅ LAMMPS .in files generated.")
PY
popd >/dev/null

echo "====================================================="
echo "  STEP 4: run LAMMPS to produce force.*"
echo "====================================================="

python3 - <<'PY'
import os, re, subprocess, shutil, sys

input_dir  = os.path.join(os.getcwd(), "lmp_format")
output_dir = os.path.join(os.getcwd(), "lmp_force_output")

lmp_exec_env = os.environ.get("LMP_EXEC", "lmp_mpi")
lmp_exec = shutil.which(lmp_exec_env) if not os.path.isabs(lmp_exec_env) else lmp_exec_env
if not lmp_exec or not os.path.exists(lmp_exec):
    print(f"❌ LAMMPS executable not found.")
    sys.exit(1)

os.makedirs(output_dir, exist_ok=True)

def extract_index(fname):
    m = re.match(r"BiSbTe_in_(\d+)\.in", fname)   # ← changed
    return int(m.group(1)) if m else -1

for f in os.listdir(output_dir):
    p = os.path.join(output_dir, f)
    if os.path.isfile(p):
        try: os.remove(p)
        except: pass

files = [f for f in os.listdir(input_dir)
         if f.startswith("BiSbTe_in_") and f.endswith(".in")]   # ← changed
files = sorted(files, key=extract_index)
if not files:
    print("❌ No BiSbTe_in_*.in found in ./lmp_format")
    sys.exit(1)

print(f"📂 Found {len(files)} jobs. Using LAMMPS: {lmp_exec}")
for fname in files:
    log_name = fname.replace(".in", ".log")
    log_path = os.path.join(output_dir, log_name)
    print(f"▶️  {fname} → {log_name}")
    cmd = [lmp_exec, "-in", fname]
    with open(log_path, "w") as logf:
        subprocess.run(cmd, stdout=logf, stderr=logf,
                       cwd=input_dir, check=False)
print("✅ LAMMPS batch finished.")
PY

echo "====================================================="
echo "  STEP 5: convert force.* → QE style qe_force.XXX.out"
echo "====================================================="

python3 - <<'PY'
import os, shutil

CONV      = 0.03889377938   # eV/Å → Ry/Bohr
ROOT      = os.getcwd()
LAMMPS_DIR = os.path.join(ROOT, "lmp_format")
STAGE_DIR  = os.path.join(ROOT, "lmp_force_output")
QE_DIR     = os.path.join(ROOT, "qe_forces")
os.makedirs(LAMMPS_DIR, exist_ok=True)
os.makedirs(STAGE_DIR,  exist_ok=True)
os.makedirs(QE_DIR,     exist_ok=True)

def list_forces(d):
    return sorted(
        [f for f in os.listdir(d)
         if f.startswith("force.") and os.path.isfile(os.path.join(d, f))],
        key=lambda s: int(s.split(".", 1)[1])
        if "." in s and s.split(".", 1)[1].isdigit() else 10**9
    )

src_files = list_forces(LAMMPS_DIR)
if not src_files:
    print("ℹ️  No force.* found yet — skip conversion.")
    raise SystemExit(0)

for f in src_files:
    shutil.copy2(os.path.join(LAMMPS_DIR, f), os.path.join(STAGE_DIR, f))
print(f"📤 Staged {len(src_files)} force.* to ./lmp_force_output")

for g in os.listdir(QE_DIR):
    p = os.path.join(QE_DIR, g)
    if os.path.isfile(p):
        try: os.remove(p)
        except: pass

def parse_force_file(path):
    dat   = []
    lines = open(path).read().splitlines()
    has_item = any(l.startswith("ITEM:") for l in lines[:20])
    start = False
    for ln in lines:
        ln = ln.strip()
        if not ln: continue
        if has_item:
            if ln.startswith("ITEM: ATOMS"):   start = True;  continue
            if ln.startswith("ITEM:") and "ATOMS" not in ln: start = False; continue
            if not start: continue
        parts = ln.split()
        if len(parts) < 8: continue
        try:
            i  = int(float(parts[0]))
            t  = int(float(parts[1]))
            fx = float(parts[5]) * CONV
            fy = float(parts[6]) * CONV
            fz = float(parts[7]) * CONV
            dat.append((i, t, fx, fy, fz))
        except: continue
    dat.sort(key=lambda x: x[0])
    return dat

def write_qe(path, dat):
    with open(path, "w") as f:
        f.write("  Forces acting on atoms (cartesian axes, Ry/au):\n\n")
        for k, (i, t, fx, fy, fz) in enumerate(dat, start=1):
            f.write(f"     atom {k:4d} type {t:2d}   "
                    f"force = {fx:12.8f} {fy:12.8f} {fz:12.8f}\n")

for f in list_forces(STAGE_DIR):
    inp  = os.path.join(STAGE_DIR, f)
    tail = f.split(".", 1)[1] if "." in f else "000"
    out  = os.path.join(QE_DIR, f"qe_force.{tail}.out")
    data = parse_force_file(inp)
    if not data:
        print(f"⚠️  {f}: no atoms parsed — skipped"); continue
    write_qe(out, data)
    print(f"✅ {f} → {os.path.basename(out)}")
PY

echo "====================================================="
echo "  STEP 6: phonopy FC2 + dispersion using band.conf"
echo "====================================================="

safe_conda_activate() {
  local had_e=0 had_u=0
  [[ $- == *e* ]] && had_e=1
  [[ $- == *u* ]] && had_u=1
  set +eu
  source "${CONDA_PATH}"
  conda deactivate >/dev/null 2>&1 || true
  conda activate "${CONDA_ENV}"
  (( had_e )) && set -e
  (( had_u )) && set -u
}

echo ">>> Activating conda env: ${CONDA_ENV}"
safe_conda_activate
command -v phonopy >/dev/null 2>&1 || { echo '❌ phonopy not found in env'; exit 1; }

[[ -d "${QE_FORCES_DIR}" ]]             || { echo "❌ ${QE_FORCES_DIR} not found"; exit 1; }
[[ -f "${INPUT_DIR}/phonopy_disp.yaml" ]] || { echo "❌ phonopy_disp.yaml not found"; exit 1; }
ls "${QE_FORCES_DIR}"/qe_force.*.out >/dev/null 2>&1 || {
  echo "❌ No qe_force.XXX.out in ${QE_FORCES_DIR}"; exit 1; }
[[ -f "${QE_CELL}" ]]  || { echo "❌ ${QE_CELL} not found"; exit 1; }
[[ -f "band.conf" ]]   || { echo "❌ band.conf not found"; exit 1; }

echo ">>> Copying needed files into ${QE_FORCES_DIR}/"
cp -f "${INPUT_DIR}/phonopy_disp.yaml" "${QE_FORCES_DIR}/"
cp -f "${QE_CELL}"  "${QE_FORCES_DIR}/"
cp -f "band.conf"   "${QE_FORCES_DIR}/"

cd "${QE_FORCES_DIR}"

echo ">>> Contents of $(pwd):"
ls

echo ">>> 6.1 → Convert qe_force.*.out → FORCE_SETS"
phonopy --qe -f qe_force.*.out

echo ">>> 6.2 → Run band structure using band.conf"
phonopy --qe -c espresso_frac.pwi -p -s band.conf

echo "🎉 Dispersion done in $(pwd)"
echo "   Check for:"
echo "     FORCE_SETS"
echo "     band.yaml"
echo "     band.pdf"
