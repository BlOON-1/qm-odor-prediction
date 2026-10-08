#!/bin/bash
#SBATCH -J ad4_hp64
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH -t 120:00:00
#SBATCH -o slurm-%x-%j.out
#SBATCH -e slurm-%x-%j.err

# AutoDock4 High-Precision Batch Pipeline.
#
# Core design:
#   1. Discover all CAS/receptor/pocket combinations automatically.
#   2. Prepare each ligand once and each complete receptor once.
#   3. Run AutoGrid once for every CAS-receptor-pocket combination.
#   4. Split exactly 100 Lamarckian-GA runs per pocket into reproducible shards.
#   5. Execute independent shards concurrently within the Slurm CPU allocation.
#   6. Resolve expected SDF/PDB filenames case-insensitively when Linux case differs.
#   7. Aggregate every DLG into machine-readable tables and a clear master report.
#
# Expected input tree:
# project_root/
# ├── autogrid4
# ├── autodock4
# ├── <CAS>/
# │   ├── <CAS>.sdf
# │   └── <receptor>/
# │       ├── <receptor>_P_0/
# │       │   ├── <receptor>.pdb
# │       │   └── <receptor>_P_0.pdb
# │       └── ...
# └── this_script.sh
#
# Recommended first check:
#   RUN_DOCKING=0 bash pipeline_scripts/run_autodock4_batch_pipeline_64core.sh
#
# Formal Slurm run:
#   sbatch pipeline_scripts/run_autodock4_batch_pipeline_64core.sh
#
# Resume the same parameter set:
#   SKIP_EXISTING=1 sbatch pipeline_scripts/run_autodock4_batch_pipeline_64core.sh

set -uo pipefail
IFS=$'\n\t'
shopt -s nullglob

# =============================================================================
# 0. Working directory and user-adjustable parameters
# =============================================================================
resolve_project_root() {
    local candidate
    for candidate in "${PIPELINE_PROJECT_DIR:-}" "${PROJECT_ROOT:-}" "${SLURM_SUBMIT_DIR:-}" "$PWD"; do
        [ -n "$candidate" ] || continue
        if [ -f "$candidate/pipeline_scripts/pipeline_config.sh" ]; then
            (cd "$candidate" && pwd -P)
            return 0
        fi
    done
    return 1
}
PROJECT_ROOT="$(resolve_project_root)" || {
    echo "ERROR: cannot resolve project root from PIPELINE_PROJECT_DIR, PROJECT_ROOT, SLURM_SUBMIT_DIR, or PWD" >&2
    exit 1
}
PIPELINE_PROJECT_DIR="$PROJECT_ROOT"
SCRIPT_DIR="$PROJECT_ROOT/pipeline_scripts"
export PROJECT_ROOT PIPELINE_PROJECT_DIR
# shellcheck source=pipeline_config.sh
source "$SCRIPT_DIR/pipeline_config.sh"
# shellcheck source=pipeline_status.sh
source "$SCRIPT_DIR/pipeline_status.sh"
WORKDIR="${WORKDIR:-$PROJECT_ROOT}"
cd "$WORKDIR" || { echo "ERROR: cannot enter WORKDIR: $WORKDIR" >&2; exit 1; }
WORKDIR="$(pwd -P)"

# Software paths are required from pipeline_config.sh or exported environment variables.
AUTOGRID="${AUTOGRID:-$WORKDIR/autogrid4}"
AUTODOCK="${AUTODOCK:-$WORKDIR/autodock4}"
PYTHON3="${PYTHON3:-python3}"  # Python 3.9+; tested logic is compatible with Python 3.9.12

# Grid-box parameters.
SPACING="${SPACING:-0.375}"
PADDING="${PADDING:-5.0}"
CENTER_METHOD="${CENTER_METHOD:-bbox}"       # bbox or centroid
ALLOW_LARGE_NPTS="${ALLOW_LARGE_NPTS:-0}"   # 0: reject npts > 126

# High-search-effort AutoDock4 parameters.
# GA_RUN_TOTAL is the exact total number of independent LGA runs per pocket.
GA_RUN_TOTAL="${GA_RUN_TOTAL:-100}"
GA_NUM_EVALS="${GA_NUM_EVALS:-10000000}"
GA_NUM_GENERATIONS="${GA_NUM_GENERATIONS:-50000}"
GA_POP_SIZE="${GA_POP_SIZE:-300}"
RMSTOL="${RMSTOL:-2.0}"

# Parallel execution.
ALLOCATED_CPUS="${SLURM_CPUS_PER_TASK:-${DOCK_CPUS:-64}}"
MAX_PARALLEL="${MAX_PARALLEL:-${DOCK_CPUS:-64}}"
if [[ ! "$ALLOCATED_CPUS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: SLURM_CPUS_PER_TASK/allocated CPU value is not a positive integer: $ALLOCATED_CPUS" >&2
    exit 1
fi
if [ "$MAX_PARALLEL" -gt "$ALLOCATED_CPUS" ]; then
    echo "WARNING: MAX_PARALLEL=$MAX_PARALLEL exceeds allocated CPUs=$ALLOCATED_CPUS; capping it." >&2
    MAX_PARALLEL="$ALLOCATED_CPUS"
fi
if [ "$MAX_PARALLEL" -gt 64 ]; then
    echo "WARNING: this workflow is capped at 64 AutoDock workers; capping MAX_PARALLEL." >&2
    MAX_PARALLEL=64
fi
# auto: choose enough shards per pocket to create about MAX_PARALLEL concurrent jobs.
# Integer override example: SHARDS_PER_POCKET=10
SHARDS_PER_POCKET="${SHARDS_PER_POCKET:-auto}"

# Execution/filtering switches.
RUN_DOCKING="${RUN_DOCKING:-1}"
SKIP_EXISTING="${SKIP_EXISTING:-1}"
OVERWRITE_PREP="${OVERWRITE_PREP:-0}"
STRICT_CAS_NAME="${STRICT_CAS_NAME:-1}"
CAS_FILTER="${CAS_FILTER:-*}"
RECEPTOR_FILTER="${RECEPTOR_FILTER:-*}"
POCKET_FILTER="${POCKET_FILTER:-*}"

# Prevent hidden multithreading inside each serial AutoDock process.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

for integer_value in "$GA_RUN_TOTAL" "$GA_NUM_EVALS" "$GA_NUM_GENERATIONS" "$GA_POP_SIZE" "$MAX_PARALLEL"; do
    if [[ ! "$integer_value" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: expected a positive integer, received: $integer_value" >&2
        exit 1
    fi
done
if [ "$SHARDS_PER_POCKET" != "auto" ] && [[ ! "$SHARDS_PER_POCKET" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: SHARDS_PER_POCKET must be 'auto' or a positive integer." >&2
    exit 1
fi
if [ "$CENTER_METHOD" != "bbox" ] && [ "$CENTER_METHOD" != "centroid" ]; then
    echo "ERROR: CENTER_METHOD must be bbox or centroid." >&2
    exit 1
fi

RESULT_ROOT="${RESULT_ROOT:-$WORKDIR/_autodock4_batch}"
RUN_TAG="${RUN_TAG:-${SLURM_JOB_ID:-local}_$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="$RESULT_ROOT/runs/$RUN_TAG"
PREP_ROOT="$RESULT_ROOT/prepared"
STATUS_ROOT="$RUN_DIR/task_status"
mkdir -p "$RUN_DIR" "$PREP_ROOT" "$STATUS_ROOT/grid" "$STATUS_ROOT/docking"

PARAM_TEXT="spacing=$SPACING|padding=$PADDING|center=$CENTER_METHOD|ga_runs=$GA_RUN_TOTAL|evals=$GA_NUM_EVALS|generations=$GA_NUM_GENERATIONS|population=$GA_POP_SIZE|rmstol=$RMSTOL"
PARAM_HASH="$(printf '%s' "$PARAM_TEXT" | sha256sum | awk '{print $1}')"
PARAM_ID="hp_${PARAM_HASH:0:12}"

PIPELINE_RUN_DIR="${PIPELINE_RUN_DIR:-$PROJECT_ROOT/_full_docking_analysis/runs/${PIPELINE_ID:-$RUN_TAG}}"
mkdir -p "$PIPELINE_RUN_DIR/logs" "$PIPELINE_RUN_DIR/stage_status"
RUN_LOG="$RUN_DIR/run.log"
ANALYSIS_DOCK_LOG="$PIPELINE_RUN_DIR/logs/docking.log"
DOCK_STAGE_START="$(date -Is)"
FAILURES="$RUN_DIR/failures.tsv"
MANIFEST="$RUN_DIR/job_manifest.tsv"
POCKET_QUEUE="$RUN_DIR/pocket_queue.tsv"
GRID_RESULTS="$RUN_DIR/grid_task_results.tsv"
SHARD_QUEUE="$RUN_DIR/docking_shard_queue.tsv"
SHARD_RESULTS="$RUN_DIR/docking_shard_results.tsv"
RECEPTOR_MD5="$RUN_DIR/receptor_md5.tsv"
PARAMETERS="$RUN_DIR/run_parameters.tsv"
CASE_MATCHES="$RUN_DIR/case_insensitive_matches.tsv"

printf 'cas\treceptor\tpocket\tstage\tmessage\n' > "$FAILURES"
printf 'queue_id\tcas\treceptor\tpocket\tpocket_pdb\tfull_receptor_pdb\tligand_sdf\tparameter_dir\tgridcenter\tnpts\tspacing\tpadding\tgpf\tdpf\tglg\tstatus\n' > "$MANIFEST"
printf 'queue_id\tcas\treceptor\tpocket\tpocket_pdb\tfull_receptor_pdb\tligand_sdf\tparameter_dir\tgpf\tdpf\tglg\n' > "$POCKET_QUEUE"
printf 'queue_id\tcas\treceptor\tpocket\tparameter_dir\tstatus\tmessage\tglg\n' > "$GRID_RESULTS"
printf 'task_id\tqueue_id\tcas\treceptor\tpocket\tparameter_dir\tshard_dir\tshard_index\tshard_count\truns_requested\tseed1\tseed2\tdpf\tdlg\tlog\n' > "$SHARD_QUEUE"
printf 'task_id\tqueue_id\tcas\treceptor\tpocket\tshard_index\truns_requested\tstatus\tmessage\tdlg\tlog\n' > "$SHARD_RESULTS"
printf 'cas\treceptor\tpocket\tmd5\tfull_receptor_pdb\n' > "$RECEPTOR_MD5"
printf 'parameter\tvalue\n' > "$PARAMETERS"
printf 'cas\treceptor\tpocket\tinput_type\texpected_path\tresolved_path\n' > "$CASE_MATCHES"

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$RUN_LOG" "$ANALYSIS_DOCK_LOG"
}

write_docking_stage_status() {
    local rc="${1:-1}" status="SUCCESS" message="Docking stage completed"
    if [ "$rc" -ne 0 ]; then
        status="FAILED"
        message="Docking stage exited with code $rc; inspect $ANALYSIS_DOCK_LOG and $FAILURES"
    fi
    pipeline_write_stage_status docking "$status" "$rc" "$message" "$WORKDIR" "$RUN_DIR" "$DOCK_STAGE_START" "$(date -Is)"
}
on_docking_exit() {
    local rc="${1:-1}"
    trap - EXIT
    write_docking_stage_status "$rc"
}
trap 'rc=$?; on_docking_exit "$rc"' EXIT
pipeline_write_stage_status docking RUNNING "" "Docking stage is running" "$WORKDIR" "$RUN_DIR" "$DOCK_STAGE_START" ""

sanitize_field() {
    printf '%s' "$1" | tr '\t\r\n' '   '
}

record_failure() {
    local cas="$1" receptor="$2" pocket="$3" stage="$4" message="$5"
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$cas" "$receptor" "$pocket" "$stage" "$(sanitize_field "$message")" >> "$FAILURES"
    log "FAIL [$cas | $receptor | $pocket | $stage] $message"
}

# Resolve one expected non-empty file in a directory.
#
# Resolution order:
#   1. Use the exact case-sensitive path when it exists and is non-empty.
#   2. Otherwise search only for the same basename ignoring letter case.
#   3. Accept exactly one non-empty case-insensitive match.
#   4. Reject zero matches or multiple matches; never guess between ambiguous files.
#
# Arguments:
#   directory expected_basename output_variable cas receptor pocket input_type
resolve_case_insensitive_file() {
    local directory="$1"
    local expected_basename="$2"
    local output_variable="$3"
    local cas="$4"
    local receptor="$5"
    local pocket="$6"
    local input_type="$7"

    local exact_path="$directory/$expected_basename"
    local -a matches=()
    local match candidate_list

    while IFS= read -r -d '' match; do
        [ -s "$match" ] && matches+=("$match")
    done < <(
        find "$directory" -maxdepth 1 -type f \
            -iname "$expected_basename" -print0 2>/dev/null | sort -z
    )

    if [ "${#matches[@]}" -eq 1 ]; then
        printf -v "$output_variable" '%s' "${matches[0]}"
        if [ "${matches[0]}" != "$exact_path" ]; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$cas" "$receptor" "$pocket" "$input_type" \
                "$(sanitize_field "$exact_path")" "$(sanitize_field "${matches[0]}")" \
                >> "$CASE_MATCHES"
            log "CASE-INSENSITIVE MATCH [$cas | $receptor | $pocket | $input_type]: expected $exact_path ; using ${matches[0]}"
        fi
        return 0
    fi

    if [ "${#matches[@]}" -eq 0 ]; then
        return 1
    fi

    candidate_list="$(printf '%s; ' "${matches[@]}")"
    log "AMBIGUOUS CASE-INSENSITIVE MATCH [$cas | $receptor | $pocket | $input_type]: expected $exact_path ; candidates: $candidate_list"
    return 2
}

is_complete_dlg() {
    local dlg="$1"
    [ -s "$dlg" ] && grep -qiE 'Successful Completion|DOCKED:.*Estimated Free Energy of Binding' "$dlg"
}

wait_for_slot() {
    while [ "$(jobs -pr | wc -l)" -ge "$MAX_PARALLEL" ]; do
        wait -n 2>/dev/null || true
    done
}

wait_for_all() {
    while [ "$(jobs -pr | wc -l)" -gt 0 ]; do
        wait -n 2>/dev/null || true
    done
}

for kv in \
    "workdir=$WORKDIR" \
    "host=$(hostname)" \
    "slurm_job_id=${SLURM_JOB_ID:-NA}" \
    "slurm_cpus_per_task=${SLURM_CPUS_PER_TASK:-NA}" \
    "max_parallel=$MAX_PARALLEL" \
    "parameter_id=$PARAM_ID" \
    "spacing_A=$SPACING" \
    "padding_each_side_A=$PADDING" \
    "center_method=$CENTER_METHOD" \
    "ga_run_total_per_pocket=$GA_RUN_TOTAL" \
    "ga_num_evals_per_run=$GA_NUM_EVALS" \
    "ga_num_generations_per_run=$GA_NUM_GENERATIONS" \
    "ga_population_size=$GA_POP_SIZE" \
    "rmstol_A=$RMSTOL" \
    "shards_per_pocket_requested=$SHARDS_PER_POCKET" \
    "run_docking=$RUN_DOCKING" \
    "skip_existing=$SKIP_EXISTING"; do
    printf '%s\t%s\n' "${kv%%=*}" "${kv#*=}" >> "$PARAMETERS"
done

# =============================================================================
# 1. Validate software
# =============================================================================
log "============================================================"
log "AutoDock4 High-Precision Batch Pipeline"
log "Work directory       : $WORKDIR"
log "Run directory        : $RUN_DIR"
log "Parameter ID         : $PARAM_ID"
log "Requested CPU workers: $MAX_PARALLEL"
log "Total LGA runs/pocket: $GA_RUN_TOTAL"
log "Evaluations/run      : $GA_NUM_EVALS"
log "Generations/run      : $GA_NUM_GENERATIONS"
log "Population size      : $GA_POP_SIZE"
log "RUN_DOCKING          : $RUN_DOCKING"
log "============================================================"

missing_software=0
for f in "$OBABEL" "$PYTHONSH" "$PREP_LIGAND" "$PREP_RECEPTOR" "$AUTOGRID" "$AUTODOCK"; do
    if [ ! -e "$f" ]; then
        log "ERROR: required software file not found: $f"
        missing_software=1
    fi
done
if ! command -v "$PYTHON3" >/dev/null 2>&1 && [ ! -x "$PYTHON3" ]; then
    log "ERROR: Python interpreter is not available: $PYTHON3"
    missing_software=1
fi
if [ "$missing_software" -ne 0 ]; then
    exit 1
fi

if ! "$PYTHON3" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 9) else 1)'; then
    log "ERROR: this script requires Python >= 3.9. Detected: $("$PYTHON3" --version 2>&1)"
    exit 1
fi
PYTHON_VERSION="$("$PYTHON3" -c 'import platform; print(platform.python_version())')"
printf 'python3_executable\t%s\n' "$PYTHON3" >> "$PARAMETERS"
printf 'python3_version\t%s\n' "$PYTHON_VERSION" >> "$PARAMETERS"
log "Python interpreter    : $PYTHON3"
log "Python version        : $PYTHON_VERSION"

chmod +x "$AUTOGRID" "$AUTODOCK" 2>/dev/null || true
"$OBABEL" -V >> "$RUN_LOG" 2>&1 || true

# =============================================================================
# 2. Ligand/receptor preparation helpers
# =============================================================================
prepare_ligand() {
    local cas="$1" ligand_sdf="$2" cas_prep="$3"
    local mol2="$cas_prep/ligand.mol2"
    local pdbqt="$cas_prep/ligand.pdbqt"
    local prep_log="$cas_prep/prepare_ligand.log"

    mkdir -p "$cas_prep"
    if [ "$OVERWRITE_PREP" != "1" ] && [ -s "$pdbqt" ]; then
        log "Reuse ligand PDBQT: $pdbqt"
        return 0
    fi

    : > "$prep_log"
    log "Prepare ligand for CAS $cas"
    if ! "$OBABEL" "$ligand_sdf" -O "$mol2" -h >> "$prep_log" 2>&1; then
        record_failure "$cas" "-" "-" "prepare_ligand_obabel" "Open Babel failed; see $prep_log"
        return 1
    fi
    if [ ! -s "$mol2" ]; then
        record_failure "$cas" "-" "-" "prepare_ligand_obabel" "MOL2 was not generated: $mol2"
        return 1
    fi
    if ! (
        cd "$cas_prep" &&
        "$PYTHONSH" "$PREP_LIGAND" -l ligand.mol2 -o ligand.pdbqt -A hydrogens
    ) >> "$prep_log" 2>&1; then
        record_failure "$cas" "-" "-" "prepare_ligand_pdbqt" "prepare_ligand4.py failed; see $prep_log"
        return 1
    fi
    if [ ! -s "$pdbqt" ]; then
        record_failure "$cas" "-" "-" "prepare_ligand_pdbqt" "Ligand PDBQT was not generated: $pdbqt"
        return 1
    fi
    return 0
}

prepare_receptor() {
    local cas="$1" receptor="$2" full_pdb="$3" receptor_prep="$4"
    local raw="$receptor_prep/${receptor}_raw.pdb"
    local clean="$receptor_prep/${receptor}_clean.pdb"
    local pdbqt="$receptor_prep/receptor.pdbqt"
    local prep_log="$receptor_prep/prepare_receptor.log"

    mkdir -p "$receptor_prep"
    if [ "$OVERWRITE_PREP" != "1" ] && [ -s "$pdbqt" ]; then
        log "Reuse receptor PDBQT: $pdbqt"
        return 0
    fi

    : > "$prep_log"
    log "Prepare receptor $receptor for CAS $cas"
    cp -f "$full_pdb" "$raw"
    grep -E '^(ATOM|TER|END)' "$raw" > "$clean" || true
    if ! grep -q '^ATOM' "$clean"; then
        record_failure "$cas" "$receptor" "-" "clean_receptor" "No ATOM records after cleaning: $clean"
        return 1
    fi
    if ! (
        cd "$receptor_prep" &&
        "$PYTHONSH" "$PREP_RECEPTOR" \
            -r "${receptor}_clean.pdb" \
            -o receptor.pdbqt \
            -A checkhydrogens \
            -U nphs_lps_waters
    ) >> "$prep_log" 2>&1; then
        record_failure "$cas" "$receptor" "-" "prepare_receptor_pdbqt" "prepare_receptor4.py failed; see $prep_log"
        return 1
    fi
    if [ ! -s "$pdbqt" ]; then
        record_failure "$cas" "$receptor" "-" "prepare_receptor_pdbqt" "Receptor PDBQT was not generated: $pdbqt"
        return 1
    fi
    return 0
}

# =============================================================================
# 3. Generate GPF/DPF and metadata for one pocket
# =============================================================================
generate_configs() {
    local cas="$1" receptor="$2" pocket="$3" pocket_pdb="$4" full_pdb="$5"
    local ligand_sdf="$6" ligand_pdbqt="$7" receptor_pdbqt="$8" outdir="$9"

    "$PYTHON3" - \
        "$cas" "$receptor" "$pocket" "$pocket_pdb" "$full_pdb" "$ligand_sdf" \
        "$ligand_pdbqt" "$receptor_pdbqt" "$outdir" \
        "$SPACING" "$PADDING" "$CENTER_METHOD" \
        "$GA_RUN_TOTAL" "$GA_NUM_EVALS" "$GA_NUM_GENERATIONS" "$GA_POP_SIZE" "$RMSTOL" \
        "$ALLOW_LARGE_NPTS" "$PARAM_ID" <<'PYEOF'
import csv
import math
import shutil
import sys
from pathlib import Path

(
    cas, receptor, pocket, pocket_pdb, full_pdb, ligand_sdf,
    ligand_pdbqt, receptor_pdbqt, outdir,
    spacing, padding, center_method,
    ga_run_total, ga_num_evals, ga_num_generations, ga_pop_size, rmstol,
    allow_large_npts, parameter_id,
) = sys.argv[1:]

spacing = float(spacing)
padding = float(padding)
ga_run_total = int(ga_run_total)
ga_num_evals = int(ga_num_evals)
ga_num_generations = int(ga_num_generations)
ga_pop_size = int(ga_pop_size)
rmstol = float(rmstol)
allow_large_npts = allow_large_npts == "1"

outdir = Path(outdir)
outdir.mkdir(parents=True, exist_ok=True)
pocket_pdb = Path(pocket_pdb)
ligand_pdbqt = Path(ligand_pdbqt)
receptor_pdbqt = Path(receptor_pdbqt)
prefix = pocket


def read_coords(path):
    coords = []
    with open(path, "r", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except ValueError:
                    continue
    if not coords:
        raise ValueError(f"No ATOM/HETATM coordinates found in {path}")
    return coords


def even_ceil(value):
    n = math.ceil(value)
    return n if n % 2 == 0 else n + 1


def unique(values):
    output, seen = [], set()
    for value in values:
        if value and value not in seen:
            output.append(value)
            seen.add(value)
    return output


def atom_types(path):
    values = []
    with open(path, "r", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                fields = line.split()
                if fields:
                    values.append(fields[-1])
    values = unique(values)
    if not values:
        raise ValueError(f"No AutoDock atom types found in {path}")
    return values


def torsdof(path):
    with open(path, "r", errors="ignore") as handle:
        for line in handle:
            if line.startswith("TORSDOF"):
                fields = line.split()
                if len(fields) >= 2:
                    return int(float(fields[1]))
    raise ValueError(f"No TORSDOF line found in ligand PDBQT: {path}")


def bbox_center(path):
    coords = read_coords(path)
    xs, ys, zs = zip(*coords)
    return ((min(xs)+max(xs))/2.0, (min(ys)+max(ys))/2.0, (min(zs)+max(zs))/2.0)

coords = read_coords(pocket_pdb)
xs, ys, zs = zip(*coords)
xmin, xmax = min(xs), max(xs)
ymin, ymax = min(ys), max(ys)
zmin, zmax = min(zs), max(zs)
centroid = (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))
bbox = ((xmin+xmax)/2.0, (ymin+ymax)/2.0, (zmin+zmax)/2.0)
center = bbox if center_method == "bbox" else centroid
lengths = (xmax-xmin, ymax-ymin, zmax-zmin)
npts = tuple(even_ceil((length + 2.0*padding)/spacing) for length in lengths)
actual_box = tuple(value*spacing for value in npts)
if max(npts) > 126 and not allow_large_npts:
    print(f"ERROR_NPTS: npts={npts[0]} {npts[1]} {npts[2]} exceeds 126", file=sys.stderr)
    sys.exit(20)

ligand_types = atom_types(ligand_pdbqt)
receptor_types = atom_types(receptor_pdbqt)
ligand_torsdof = torsdof(ligand_pdbqt)
about = bbox_center(ligand_pdbqt)

shutil.copy2(ligand_pdbqt, outdir / "ligand.pdbqt")
shutil.copy2(receptor_pdbqt, outdir / "receptor.pdbqt")

gpf = outdir / f"{prefix}.gpf"
dpf = outdir / f"{prefix}.base.dpf"
glg = outdir / f"{prefix}.glg"
box_tsv = outdir / f"{prefix}_box.tsv"
meta_tsv = outdir / f"{prefix}_config.tsv"

with open(gpf, "w") as handle:
    handle.write(f"# CAS={cas}; receptor={receptor}; pocket={pocket}; parameter_id={parameter_id}\n")
    handle.write(f"npts {npts[0]} {npts[1]} {npts[2]}\n")
    handle.write(f"gridfld {prefix}.maps.fld\n")
    handle.write(f"spacing {spacing:.3f}\n")
    handle.write("receptor_types " + " ".join(receptor_types) + "\n")
    handle.write("ligand_types " + " ".join(ligand_types) + "\n")
    handle.write("receptor receptor.pdbqt\n")
    handle.write(f"gridcenter {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n")
    handle.write("smooth 0.5\n")
    for atom_type in ligand_types:
        handle.write(f"map {prefix}.{atom_type}.map\n")
    handle.write(f"elecmap {prefix}.e.map\n")
    handle.write(f"dsolvmap {prefix}.d.map\n")
    handle.write("dielectric -0.1465\n")

with open(dpf, "w") as handle:
    handle.write(f"# Base DPF; shards replace seed and ga_run. CAS={cas}; receptor={receptor}; pocket={pocket}\n")
    handle.write("autodock_parameter_version 4.2\n")
    handle.write("outlev 1\n")
    handle.write("intelec\n")
    handle.write("seed 1 1\n")
    handle.write("ligand_types " + " ".join(ligand_types) + "\n")
    handle.write(f"fld {prefix}.maps.fld\n")
    for atom_type in ligand_types:
        handle.write(f"map {prefix}.{atom_type}.map\n")
    handle.write(f"elecmap {prefix}.e.map\n")
    handle.write(f"desolvmap {prefix}.d.map\n")
    handle.write("move ligand.pdbqt\n")
    handle.write(f"about {about[0]:.3f} {about[1]:.3f} {about[2]:.3f}\n")
    handle.write("tran0 random\n")
    handle.write("quaternion0 random\n")
    handle.write("dihe0 random\n")
    handle.write(f"torsdof {ligand_torsdof}\n")
    handle.write(f"rmstol {rmstol:.1f}\n")
    handle.write("extnrg 1000.0\n")
    handle.write("e0max 0.0 10000\n")
    handle.write(f"ga_pop_size {ga_pop_size}\n")
    handle.write(f"ga_num_evals {ga_num_evals}\n")
    handle.write(f"ga_num_generations {ga_num_generations}\n")
    handle.write("ga_elitism 1\n")
    handle.write("ga_mutation_rate 0.02\n")
    handle.write("ga_crossover_rate 0.8\n")
    handle.write("ga_window_size 10\n")
    handle.write("ga_cauchy_alpha 0.0\n")
    handle.write("ga_cauchy_beta 1.0\n")
    handle.write("set_ga\n")
    handle.write("sw_max_its 300\n")
    handle.write("sw_max_succ 4\n")
    handle.write("sw_max_fail 4\n")
    handle.write("sw_rho 1.0\n")
    handle.write("sw_lb_rho 0.01\n")
    handle.write("ls_search_freq 0.06\n")
    handle.write("set_psw1\n")
    handle.write("unbound_model bound\n")
    handle.write(f"ga_run {ga_run_total}\n")
    handle.write("analysis\n")

box_fields = [
    "cas", "receptor", "pocket", "pocket_pdb", "atom_count", "center_method",
    "center_x", "center_y", "center_z", "xmin", "xmax", "ymin", "ymax",
    "zmin", "zmax", "spacing", "padding_each_side", "npts_x", "npts_y",
    "npts_z", "actual_box_x", "actual_box_y", "actual_box_z",
]
box_row = {
    "cas": cas, "receptor": receptor, "pocket": pocket, "pocket_pdb": str(pocket_pdb),
    "atom_count": len(coords), "center_method": center_method,
    "center_x": f"{center[0]:.3f}", "center_y": f"{center[1]:.3f}", "center_z": f"{center[2]:.3f}",
    "xmin": f"{xmin:.3f}", "xmax": f"{xmax:.3f}", "ymin": f"{ymin:.3f}",
    "ymax": f"{ymax:.3f}", "zmin": f"{zmin:.3f}", "zmax": f"{zmax:.3f}",
    "spacing": f"{spacing:.3f}", "padding_each_side": f"{padding:.3f}",
    "npts_x": npts[0], "npts_y": npts[1], "npts_z": npts[2],
    "actual_box_x": f"{actual_box[0]:.3f}", "actual_box_y": f"{actual_box[1]:.3f}",
    "actual_box_z": f"{actual_box[2]:.3f}",
}
with open(box_tsv, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=box_fields, delimiter="\t")
    writer.writeheader(); writer.writerow(box_row)

meta_fields = [
    "cas", "receptor", "pocket", "pocket_pdb", "full_receptor_pdb", "ligand_sdf",
    "parameter_dir", "gridcenter", "npts", "spacing", "padding", "gpf", "dpf", "glg", "status",
]
meta_row = {
    "cas": cas, "receptor": receptor, "pocket": pocket, "pocket_pdb": str(pocket_pdb),
    "full_receptor_pdb": full_pdb, "ligand_sdf": ligand_sdf,
    "parameter_dir": str(outdir),
    "gridcenter": f"{center[0]:.3f} {center[1]:.3f} {center[2]:.3f}",
    "npts": f"{npts[0]} {npts[1]} {npts[2]}", "spacing": f"{spacing:.3f}",
    "padding": f"{padding:.3f}", "gpf": str(gpf), "dpf": str(dpf), "glg": str(glg),
    "status": "CONFIGURED",
}
with open(meta_tsv, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=meta_fields, delimiter="\t")
    writer.writeheader(); writer.writerow(meta_row)
print(meta_tsv)
PYEOF
}

# =============================================================================
# 4. Discover all inputs and build the pocket queue
# =============================================================================
cas_dirs=("$WORKDIR"/*-*-*)
if [ "${#cas_dirs[@]}" -eq 0 ]; then
    log "ERROR: no immediate subdirectories matching *-*-* were found under $WORKDIR"
    exit 1
fi

cas_seen=0
receptor_seen=0
pocket_seen=0
configured_count=0
failed_count=0
queue_index=0

for cas_dir in "${cas_dirs[@]}"; do
    [ -d "$cas_dir" ] || continue
    cas="$(basename "$cas_dir")"
    if [ "$STRICT_CAS_NAME" = "1" ] && [[ ! "$cas" =~ ^[0-9]+-[0-9]+-[0-9]+$ ]]; then
        continue
    fi
    [[ "$cas" == $CAS_FILTER ]] || continue
    cas_seen=$((cas_seen + 1))

    ligand_sdf=""
    resolve_case_insensitive_file         "$cas_dir" "$cas.sdf" ligand_sdf         "$cas" "-" "-" "ligand_sdf"
    ligand_rc=$?
    if [ "$ligand_rc" -ne 0 ]; then
        if [ "$ligand_rc" -eq 2 ]; then
            message="Ambiguous case-insensitive ligand SDF matches for expected path: $cas_dir/$cas.sdf"
        else
            message="Missing or empty ligand SDF, including case-insensitive check: $cas_dir/$cas.sdf"
        fi
        record_failure "$cas" "-" "-" "input_discovery" "$message"
        failed_count=$((failed_count + 1))
        continue
    fi

    cas_prep="$PREP_ROOT/$cas"
    if ! prepare_ligand "$cas" "$ligand_sdf" "$cas_prep"; then
        failed_count=$((failed_count + 1))
        continue
    fi
    ligand_pdbqt="$cas_prep/ligand.pdbqt"

    receptor_dirs=("$cas_dir"/*/)
    declare -A receptor_case_counts=()
    for receptor_dir_raw in "${receptor_dirs[@]}"; do
        receptor_dir="${receptor_dir_raw%/}"
        [ -d "$receptor_dir" ] || continue
        receptor_key="${receptor_dir##*/}"
        receptor_key="${receptor_key,,}"
        receptor_case_counts["$receptor_key"]=$(( ${receptor_case_counts["$receptor_key"]:-0} + 1 ))
    done
    for receptor_dir_raw in "${receptor_dirs[@]}"; do
        receptor_dir="${receptor_dir_raw%/}"
        [ -d "$receptor_dir" ] || continue
        receptor="$(basename "$receptor_dir")"
        receptor_key="${receptor,,}"
        if [ "${receptor_case_counts["$receptor_key"]:-0}" -gt 1 ]; then
            record_failure "$cas" "$receptor" "-" "input_discovery" "Ambiguous receptor directories differing only by case under $cas_dir"
            failed_count=$((failed_count + 1))
            continue
        fi
        [[ "$receptor" == $RECEPTOR_FILTER ]] || continue
        receptor_seen=$((receptor_seen + 1))

        # Linux filenames are case-sensitive. Discover pocket directories using
        # an exact receptor_P_* pattern while ignoring letter case.
        valid_pocket_dirs=()
        while IFS= read -r -d '' pocket_dir; do
            pocket="$(basename "$pocket_dir")"
            [[ "$pocket" == $POCKET_FILTER ]] || continue
            valid_pocket_dirs+=("$pocket_dir")
            if [[ "$pocket" != "${receptor}_P_"* ]]; then
                printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
                    "$cas" "$receptor" "$pocket" "pocket_directory" \
                    "$(sanitize_field "$receptor_dir/${receptor}_P_*")" "$(sanitize_field "$pocket_dir")" \
                    >> "$CASE_MATCHES"
            fi
        done < <(
            find "$receptor_dir" -mindepth 1 -maxdepth 1 -type d                 -iname "${receptor}_P_*" -print0 2>/dev/null | sort -z
        )
        if [ "${#valid_pocket_dirs[@]}" -eq 0 ]; then
            record_failure "$cas" "$receptor" "-" "input_discovery" "No pocket folders matching ${receptor}_P_*"
            failed_count=$((failed_count + 1))
            continue
        fi
        declare -A pocket_case_counts=()
        for pocket_dir in "${valid_pocket_dirs[@]}"; do
            pocket_key="$(basename "$pocket_dir")"
            pocket_key="${pocket_key,,}"
            pocket_case_counts["$pocket_key"]=$(( ${pocket_case_counts["$pocket_key"]:-0} + 1 ))
        done
        ambiguous_pocket_dirs=0
        for pocket_key in "${!pocket_case_counts[@]}"; do
            if [ "${pocket_case_counts["$pocket_key"]}" -gt 1 ]; then
                record_failure "$cas" "$receptor" "$pocket_key" "input_discovery" "Ambiguous pocket directories differing only by case under $receptor_dir"
                ambiguous_pocket_dirs=1
            fi
        done
        if [ "$ambiguous_pocket_dirs" -ne 0 ]; then
            failed_count=$((failed_count + 1))
            continue
        fi

        full_pdb=""
        for pocket_dir in "${valid_pocket_dirs[@]}"; do
            pocket="$(basename "$pocket_dir")"
            candidate_full=""
            resolve_case_insensitive_file                 "$pocket_dir" "$receptor.pdb" candidate_full                 "$cas" "$receptor" "$pocket" "complete_receptor_pdb"
            full_rc=$?

            if [ "$full_rc" -eq 0 ]; then
                current_md5="$(md5sum "$candidate_full" | awk '{print $1}')"
                printf '%s\t%s\t%s\t%s\t%s\n'                     "$cas" "$receptor" "$pocket" "$current_md5" "$candidate_full"                     >> "$RECEPTOR_MD5"
                [ -n "$full_pdb" ] || full_pdb="$candidate_full"
            elif [ "$full_rc" -eq 2 ]; then
                record_failure "$cas" "$receptor" "$pocket" "input_discovery"                     "Ambiguous case-insensitive complete receptor PDB matches for expected path: $pocket_dir/$receptor.pdb"
            else
                record_failure "$cas" "$receptor" "$pocket" "input_discovery"                     "Missing or empty complete receptor PDB, including case-insensitive check: $pocket_dir/$receptor.pdb"
            fi
        done
        if [ -z "$full_pdb" ]; then
            failed_count=$((failed_count + 1))
            continue
        fi

        md5_unique_count="$(awk -F '\t' -v c="$cas" -v r="$receptor" 'NR>1 && $1==c && $2==r {print $4}' "$RECEPTOR_MD5" | sort -u | wc -l)"
        if [ "$md5_unique_count" -gt 1 ]; then
            log "WARNING [$cas | $receptor]: complete receptor PDB copies differ; using $full_pdb"
        fi

        receptor_prep="$PREP_ROOT/$cas/$receptor"
        if ! prepare_receptor "$cas" "$receptor" "$full_pdb" "$receptor_prep"; then
            failed_count=$((failed_count + 1))
            continue
        fi
        receptor_pdbqt="$receptor_prep/receptor.pdbqt"

        for pocket_dir in "${valid_pocket_dirs[@]}"; do
            pocket="$(basename "$pocket_dir")"
            pocket_pdb=""
            pocket_seen=$((pocket_seen + 1))

            resolve_case_insensitive_file                 "$pocket_dir" "$pocket.pdb" pocket_pdb                 "$cas" "$receptor" "$pocket" "pocket_pdb"
            pocket_rc=$?
            if [ "$pocket_rc" -ne 0 ]; then
                if [ "$pocket_rc" -eq 2 ]; then
                    message="Ambiguous case-insensitive pocket PDB matches for expected path: $pocket_dir/$pocket.pdb"
                else
                    message="Missing or empty pocket PDB, including case-insensitive check: $pocket_dir/$pocket.pdb"
                fi
                record_failure "$cas" "$receptor" "$pocket" "input_discovery" "$message"
                failed_count=$((failed_count + 1))
                continue
            fi

            queue_index=$((queue_index + 1))
            queue_id="$(printf 'Q%04d' "$queue_index")"
            parameter_dir="$pocket_dir/ad4_runs/$PARAM_ID"
            mkdir -p "$parameter_dir"
            config_stdout="$parameter_dir/${pocket}_config.stdout"
            config_stderr="$parameter_dir/${pocket}_config.stderr"

            log "Configure $queue_id: $cas | $receptor | $pocket"
            generate_configs \
                "$cas" "$receptor" "$pocket" "$pocket_pdb" "$full_pdb" "$ligand_sdf" \
                "$ligand_pdbqt" "$receptor_pdbqt" "$parameter_dir" \
                > "$config_stdout" 2> "$config_stderr"
            rc=$?
            if [ "$rc" -ne 0 ]; then
                if [ "$rc" -eq 20 ]; then
                    message="Grid npts exceeds 126; refine the pocket or adjust SPACING/PADDING. See $config_stderr"
                else
                    message="Configuration generation failed with exit code $rc. See $config_stderr"
                fi
                record_failure "$cas" "$receptor" "$pocket" "generate_configs" "$message"
                failed_count=$((failed_count + 1))
                continue
            fi

            meta_tsv="$parameter_dir/${pocket}_config.tsv"
            if [ ! -s "$meta_tsv" ]; then
                record_failure "$cas" "$receptor" "$pocket" "generate_configs" "Missing metadata: $meta_tsv"
                failed_count=$((failed_count + 1))
                continue
            fi

            IFS=$'\t' read -r _mcas _mrec _mpocket _mpocketpdb _mfull _mligand _mparam gridcenter npts mspacing mpadding gpf dpf glg mstatus \
                < <(tail -n 1 "$meta_tsv")
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$queue_id" "$cas" "$receptor" "$pocket" "$pocket_pdb" "$full_pdb" "$ligand_sdf" \
                "$parameter_dir" "$gridcenter" "$npts" "$mspacing" "$mpadding" "$gpf" "$dpf" "$glg" "CONFIGURED" >> "$MANIFEST"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$queue_id" "$cas" "$receptor" "$pocket" "$pocket_pdb" "$full_pdb" "$ligand_sdf" \
                "$parameter_dir" "$gpf" "$dpf" "$glg" >> "$POCKET_QUEUE"
            configured_count=$((configured_count + 1))
        done
    done
done

if [ "$configured_count" -eq 0 ]; then
    log "ERROR: no CAS-receptor-pocket combination was configured."
    exit 2
fi

log "Discovery complete: CAS=$cas_seen receptors=$receptor_seen pockets=$pocket_seen configured=$configured_count"

# =============================================================================
# 5. AutoGrid phase, parallel across pockets
# =============================================================================
verify_grid_outputs() {
    local gpf="$1" parameter_dir="$2"
    local output_name
    while read -r output_name; do
        [ -n "$output_name" ] || continue
        if [ ! -s "$parameter_dir/$output_name" ]; then
            return 1
        fi
    done < <(awk '$1=="map" || $1=="elecmap" || $1=="dsolvmap" || $1=="gridfld" {print $2}' "$gpf")
    return 0
}

grid_task() {
    local queue_id="$1" cas="$2" receptor="$3" pocket="$4" parameter_dir="$5" gpf="$6" glg="$7"
    local status_file="$STATUS_ROOT/grid/${queue_id}.tsv"
    local task_log="$parameter_dir/${pocket}_autogrid.pipeline.log"
    local status message

    if [ "$SKIP_EXISTING" = "1" ] && [ -s "$glg" ] && verify_grid_outputs "$gpf" "$parameter_dir"; then
        status="SKIPPED_COMPLETE"
        message="Existing complete grid outputs reused"
    else
        : > "$task_log"
        if (cd "$parameter_dir" && "$AUTOGRID" -p "$(basename "$gpf")" -l "$(basename "$glg")") >> "$task_log" 2>&1 \
            && [ -s "$glg" ] && verify_grid_outputs "$gpf" "$parameter_dir"; then
            status="SUCCESS"
            message="AutoGrid4 completed"
        else
            status="FAILED"
            message="AutoGrid4 failed or required map files are missing; inspect $task_log and $glg"
        fi
    fi
    printf 'queue_id\tcas\treceptor\tpocket\tparameter_dir\tstatus\tmessage\tglg\n' > "$status_file"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$queue_id" "$cas" "$receptor" "$pocket" "$parameter_dir" "$status" "$(sanitize_field "$message")" "$glg" >> "$status_file"
}

if [ "$RUN_DOCKING" = "1" ]; then
    log "Starting AutoGrid4 phase with up to $MAX_PARALLEL concurrent tasks"
    while IFS=$'\t' read -r queue_id cas receptor pocket pocket_pdb full_pdb ligand_sdf parameter_dir gpf dpf glg; do
        [ "$queue_id" = "queue_id" ] && continue
        wait_for_slot
        grid_task "$queue_id" "$cas" "$receptor" "$pocket" "$parameter_dir" "$gpf" "$glg" &
    done < "$POCKET_QUEUE"
    wait_for_all

    for status_file in "$STATUS_ROOT"/grid/*.tsv; do
        [ -s "$status_file" ] || continue
        tail -n 1 "$status_file" >> "$GRID_RESULTS"
    done
else
    log "RUN_DOCKING=$RUN_DOCKING: AutoGrid4 and AutoDock4 execution skipped"
fi

# =============================================================================
# 6. Build reproducible docking shards and run them within the CPU allocation
# =============================================================================
if [ "$RUN_DOCKING" = "1" ]; then
    grid_ok_count="$(awk -F '\t' 'NR>1 && ($6=="SUCCESS" || $6=="SKIPPED_COMPLETE") {n++} END{print n+0}' "$GRID_RESULTS")"
    if [ "$grid_ok_count" -eq 0 ]; then
        while IFS=$'\t' read -r queue_id cas receptor pocket parameter_dir status message glg; do
            [ "$queue_id" = "queue_id" ] && continue
            record_failure "$cas" "$receptor" "$pocket" "autogrid4" "$message"
        done < "$GRID_RESULTS"
        log "ERROR: no grid task completed successfully."
        exit 3
    fi

    if [ "$SHARDS_PER_POCKET" = "auto" ]; then
        actual_shards=$(( (MAX_PARALLEL + grid_ok_count - 1) / grid_ok_count ))
        [ "$actual_shards" -le "$GA_RUN_TOTAL" ] || actual_shards="$GA_RUN_TOTAL"
        [ "$actual_shards" -ge 1 ] || actual_shards=1
    else
        actual_shards="$SHARDS_PER_POCKET"
        [ "$actual_shards" -le "$GA_RUN_TOTAL" ] || actual_shards="$GA_RUN_TOTAL"
    fi

    printf 'actual_shards_per_pocket\t%s\n' "$actual_shards" >> "$PARAMETERS"
    printf 'planned_concurrent_docking_tasks\t%s\n' "$((grid_ok_count * actual_shards))" >> "$PARAMETERS"
    log "Grid-success pockets: $grid_ok_count"
    log "Docking shards/pocket: $actual_shards"
    log "Planned shard tasks   : $((grid_ok_count * actual_shards))"

    shard_base_runs=$((GA_RUN_TOTAL / actual_shards))
    shard_remainder=$((GA_RUN_TOTAL % actual_shards))
    task_index=0

    while IFS=$'\t' read -r queue_id cas receptor pocket parameter_dir grid_status grid_message glg; do
        [ "$queue_id" = "queue_id" ] && continue
        if [ "$grid_status" != "SUCCESS" ] && [ "$grid_status" != "SKIPPED_COMPLETE" ]; then
            record_failure "$cas" "$receptor" "$pocket" "autogrid4" "$grid_message"
            continue
        fi

        base_dpf="$parameter_dir/${pocket}.base.dpf"
        shard_dir="$parameter_dir/shards_${actual_shards}"
        mkdir -p "$shard_dir"
        ln -sfn ../ligand.pdbqt "$shard_dir/ligand.pdbqt"
        for input_file in "$parameter_dir"/"$pocket".*.map "$parameter_dir/$pocket.maps.fld"; do
            [ -e "$input_file" ] || continue
            ln -sfn "../$(basename "$input_file")" "$shard_dir/$(basename "$input_file")"
        done

        for ((shard_index=1; shard_index<=actual_shards; shard_index++)); do
            runs_requested="$shard_base_runs"
            if [ "$shard_index" -le "$shard_remainder" ]; then
                runs_requested=$((runs_requested + 1))
            fi
            task_index=$((task_index + 1))
            task_id="$(printf 'D%06d' "$task_index")"
            shard_name="$(printf 'shard_%03d' "$shard_index")"
            shard_dpf="$shard_dir/${shard_name}.dpf"
            shard_dlg="$shard_dir/${shard_name}.dlg"
            shard_log="$shard_dir/${shard_name}.pipeline.log"

            seed_hash="$(printf '%s' "$PARAM_ID|$cas|$receptor|$pocket|$actual_shards|$shard_index" | sha256sum | awk '{print $1}')"
            seed1=$((16#${seed_hash:0:7} + 1))
            seed2=$((16#${seed_hash:7:7} + 1))

            awk -v s1="$seed1" -v s2="$seed2" -v runs="$runs_requested" '
                $1=="seed"   {print "seed " s1 " " s2; next}
                $1=="ga_run" {print "ga_run " runs; next}
                {print}
            ' "$base_dpf" > "$shard_dpf"

            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$task_id" "$queue_id" "$cas" "$receptor" "$pocket" "$parameter_dir" "$shard_dir" \
                "$shard_index" "$actual_shards" "$runs_requested" "$seed1" "$seed2" "$shard_dpf" "$shard_dlg" "$shard_log" >> "$SHARD_QUEUE"
        done
    done < "$GRID_RESULTS"

    docking_task() {
        local task_id="$1" queue_id="$2" cas="$3" receptor="$4" pocket="$5"
        local shard_dir="$6" shard_index="$7" runs_requested="$8" dpf="$9" dlg="${10}" task_log="${11}"
        local status_file="$STATUS_ROOT/docking/${task_id}.tsv"
        local status message

        if [ "$SKIP_EXISTING" = "1" ] && is_complete_dlg "$dlg"; then
            status="SKIPPED_COMPLETE"
            message="Existing complete DLG reused"
        else
            : > "$task_log"
            if (cd "$shard_dir" && "$AUTODOCK" -p "$(basename "$dpf")" -l "$(basename "$dlg")") >> "$task_log" 2>&1 \
                && is_complete_dlg "$dlg"; then
                status="SUCCESS"
                message="AutoDock4 shard completed"
            else
                status="FAILED"
                message="AutoDock4 failed or DLG lacks a completion/energy marker; inspect $task_log and $dlg"
            fi
        fi

        printf 'task_id\tqueue_id\tcas\treceptor\tpocket\tshard_index\truns_requested\tstatus\tmessage\tdlg\tlog\n' > "$status_file"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$task_id" "$queue_id" "$cas" "$receptor" "$pocket" "$shard_index" "$runs_requested" \
            "$status" "$(sanitize_field "$message")" "$dlg" "$task_log" >> "$status_file"
    }

    log "Starting AutoDock4 phase with up to $MAX_PARALLEL concurrent serial processes"
    while IFS=$'\t' read -r task_id queue_id cas receptor pocket parameter_dir shard_dir shard_index shard_count runs_requested seed1 seed2 dpf dlg task_log; do
        [ "$task_id" = "task_id" ] && continue
        wait_for_slot
        docking_task "$task_id" "$queue_id" "$cas" "$receptor" "$pocket" "$shard_dir" "$shard_index" "$runs_requested" "$dpf" "$dlg" "$task_log" &
    done < "$SHARD_QUEUE"
    wait_for_all

    for status_file in "$STATUS_ROOT"/docking/*.tsv; do
        [ -s "$status_file" ] || continue
        tail -n 1 "$status_file" >> "$SHARD_RESULTS"
    done

    while IFS=$'\t' read -r task_id queue_id cas receptor pocket shard_index runs_requested status message dlg task_log; do
        [ "$task_id" = "task_id" ] && continue
        if [ "$status" != "SUCCESS" ] && [ "$status" != "SKIPPED_COMPLETE" ]; then
            record_failure "$cas" "$receptor" "$pocket" "autodock4_shard_$shard_index" "$message"
        fi
    done < "$SHARD_RESULTS"
fi

# =============================================================================
# 7. Aggregate all energies and create the master reports
# =============================================================================
"$PYTHON3" - \
    "$PARAMETERS" "$MANIFEST" "$GRID_RESULTS" "$SHARD_QUEUE" "$SHARD_RESULTS" "$FAILURES" "$CASE_MATCHES" \
    "$RUN_DIR" "$RUN_DOCKING" "$GA_RUN_TOTAL" "$PARAM_ID" <<'PYEOF'
import csv
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

(
    parameters_file, manifest_file, grid_results_file, shard_queue_file,
    shard_results_file, failures_file, case_matches_file, run_dir, run_docking,
    ga_run_total, parameter_id,
) = sys.argv[1:]

run_dir = Path(run_dir)
run_docking = run_docking == "1"
ga_run_total = int(ga_run_total)

energy_re = re.compile(
    r"Estimated Free Energy of Binding\s*=\s*([-+]?\d+(?:\.\d+)?)\s*kcal/mol",
    re.IGNORECASE,
)


def read_tsv(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def quantile(values, q):
    if not values:
        return math.nan
    values = sorted(values)
    if len(values) == 1:
        return values[0]
    pos = (len(values) - 1) * q
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return values[lo]
    return values[lo] + (values[hi] - values[lo]) * (pos - lo)


def arithmetic_mean(values):
    return sum(values) / len(values)


def median_value(values):
    ordered = sorted(values)
    n = len(ordered)
    middle = n // 2
    if n % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2.0


def sample_standard_deviation(values):
    if len(values) < 2:
        return 0.0
    mean_value = arithmetic_mean(values)
    squared_deviation_sum = sum((value - mean_value) ** 2 for value in values)
    return math.sqrt(squared_deviation_sum / (len(values) - 1))

parameters = {row["parameter"]: row["value"] for row in read_tsv(parameters_file)}
manifest_rows = read_tsv(manifest_file)
grid_rows = read_tsv(grid_results_file)
shard_queue = read_tsv(shard_queue_file)
shard_results = {row["task_id"]: row for row in read_tsv(shard_results_file)}
failure_rows = read_tsv(failures_file)
case_match_rows = read_tsv(case_matches_file)

manifest_by_key = {(r["cas"], r["receptor"], r["pocket"]): r for r in manifest_rows}
shards_by_key = defaultdict(list)
for row in shard_queue:
    shards_by_key[(row["cas"], row["receptor"], row["pocket"])].append(row)

grid_by_key = {(r["cas"], r["receptor"], r["pocket"]): r for r in grid_rows}
all_energy_rows = []
pocket_rows = []

for key, manifest in sorted(manifest_by_key.items()):
    cas, receptor, pocket = key
    shards = sorted(shards_by_key.get(key, []), key=lambda r: int(r["shard_index"]))
    requested_runs = sum(int(r["runs_requested"]) for r in shards)
    energies = []
    completed_shards = 0
    failed_shards = 0
    best_dlg = ""
    best_energy = math.inf

    for shard in shards:
        task_status = shard_results.get(shard["task_id"], {})
        status = task_status.get("status", "NOT_RUN")
        if status in {"SUCCESS", "SKIPPED_COMPLETE"}:
            completed_shards += 1
        elif status not in {"NOT_RUN", ""}:
            failed_shards += 1

        dlg = Path(shard["dlg"])
        shard_energies = []
        completion_marker = 0
        if dlg.is_file() and dlg.stat().st_size > 0:
            text = dlg.read_text(errors="ignore")
            shard_energies = [float(x) for x in energy_re.findall(text)]
            completion_marker = int("Successful Completion" in text)

        for local_index, energy in enumerate(shard_energies, start=1):
            energies.append(energy)
            all_energy_rows.append({
                "cas": cas,
                "receptor": receptor,
                "pocket": pocket,
                "shard_index": shard["shard_index"],
                "runs_requested_in_shard": shard["runs_requested"],
                "energy_record_index_in_shard": local_index,
                "estimated_binding_energy_kcal_mol": f"{energy:.3f}",
                "seed1": shard["seed1"],
                "seed2": shard["seed2"],
                "successful_completion_marker": completion_marker,
                "dlg": str(dlg),
            })
            if energy < best_energy:
                best_energy = energy
                best_dlg = str(dlg)

    grid_status = grid_by_key.get(key, {}).get("status", "NOT_RUN")
    if not run_docking:
        overall_status = "CONFIGURED_NOT_RUN"
    elif requested_runs == 0:
        overall_status = "NO_DOCKING_SHARDS"
    elif len(energies) >= requested_runs and failed_shards == 0:
        overall_status = "COMPLETE"
    elif energies:
        overall_status = "PARTIAL"
    else:
        overall_status = "FAILED"

    if energies:
        mean = arithmetic_mean(energies)
        median = median_value(energies)
        sd = sample_standard_deviation(energies)
        q25 = quantile(energies, 0.25)
        q75 = quantile(energies, 0.75)
        worst = max(energies)
    else:
        mean = median = sd = q25 = q75 = worst = math.nan

    completeness = 0.0 if requested_runs == 0 else min(100.0, 100.0 * len(energies) / requested_runs)
    pocket_rows.append({
        "cas": cas,
        "receptor": receptor,
        "pocket": pocket,
        "status": overall_status,
        "grid_status": grid_status,
        "requested_lga_runs": requested_runs if requested_runs else ga_run_total,
        "energy_records_found": len(energies),
        "completion_percent": f"{completeness:.1f}",
        "shards_expected": len(shards),
        "shards_completed_or_reused": completed_shards,
        "shards_failed": failed_shards,
        "best_energy_kcal_mol": "" if not energies else f"{best_energy:.3f}",
        "median_energy_kcal_mol": "" if not energies else f"{median:.3f}",
        "mean_energy_kcal_mol": "" if not energies else f"{mean:.3f}",
        "sd_energy_kcal_mol": "" if not energies else f"{sd:.3f}",
        "q25_energy_kcal_mol": "" if not energies else f"{q25:.3f}",
        "q75_energy_kcal_mol": "" if not energies else f"{q75:.3f}",
        "worst_energy_kcal_mol": "" if not energies else f"{worst:.3f}",
        "gridcenter": manifest["gridcenter"],
        "npts": manifest["npts"],
        "spacing_A": manifest["spacing"],
        "padding_each_side_A": manifest["padding"],
        "parameter_dir": manifest["parameter_dir"],
        "best_dlg": best_dlg,
    })

    if best_dlg:
        link = Path(manifest["parameter_dir"]) / "BEST_DLG.dlg"
        try:
            if link.exists() or link.is_symlink():
                link.unlink()
            link.symlink_to(Path(best_dlg).resolve())
        except OSError:
            pass

all_energies_path = run_dir / "all_binding_energies.tsv"
pocket_summary_path = run_dir / "pocket_summary.tsv"
top_hits_path = run_dir / "top_hits.tsv"
report_md_path = run_dir / "MASTER_REPORT.md"
report_txt_path = run_dir / "MASTER_REPORT.txt"
readme_path = run_dir / "RESULT_FILES_README.txt"

energy_fields = [
    "cas", "receptor", "pocket", "shard_index", "runs_requested_in_shard",
    "energy_record_index_in_shard", "estimated_binding_energy_kcal_mol",
    "seed1", "seed2", "successful_completion_marker", "dlg",
]
with open(all_energies_path, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=energy_fields, delimiter="\t")
    writer.writeheader(); writer.writerows(all_energy_rows)

summary_fields = [
    "cas", "receptor", "pocket", "status", "grid_status", "requested_lga_runs",
    "energy_records_found", "completion_percent", "shards_expected",
    "shards_completed_or_reused", "shards_failed", "best_energy_kcal_mol",
    "median_energy_kcal_mol", "mean_energy_kcal_mol", "sd_energy_kcal_mol",
    "q25_energy_kcal_mol", "q75_energy_kcal_mol", "worst_energy_kcal_mol",
    "gridcenter", "npts", "spacing_A", "padding_each_side_A", "parameter_dir", "best_dlg",
]
with open(pocket_summary_path, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
    writer.writeheader(); writer.writerows(pocket_rows)

ranked = sorted(
    [r for r in pocket_rows if r["best_energy_kcal_mol"]],
    key=lambda r: float(r["best_energy_kcal_mol"]),
)
with open(top_hits_path, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
    writer.writeheader(); writer.writerows(ranked)

status_counts = defaultdict(int)
for row in pocket_rows:
    status_counts[row["status"]] += 1

total_requested = sum(int(r["requested_lga_runs"]) for r in pocket_rows) if run_docking else 0
total_energy_records = sum(int(r["energy_records_found"]) for r in pocket_rows)
total_shards = sum(int(r["shards_expected"]) for r in pocket_rows)
total_failed_shards = sum(int(r["shards_failed"]) for r in pocket_rows)

lines = []
lines.append("# AutoDock4 High-Precision Batch Pipeline — Master Report")
lines.append("")
lines.append("## 1. Run overview")
lines.append("")
lines.append(f"- Parameter ID: `{parameter_id}`")
lines.append(f"- Working directory: `{parameters.get('workdir', '')}`")
lines.append(f"- Host: `{parameters.get('host', '')}`")
lines.append(f"- Slurm job ID: `{parameters.get('slurm_job_id', '')}`")
lines.append(f"- Maximum concurrent AutoDock processes: **{parameters.get('max_parallel', '')}**")
lines.append(f"- Configured CAS-receptor-pocket combinations: **{len(pocket_rows)}**")
lines.append(f"- Requested LGA runs per pocket: **{parameters.get('ga_run_total_per_pocket', ga_run_total)}**")
lines.append(f"- Total requested LGA runs: **{total_requested}**")
lines.append(f"- Energy records recovered: **{total_energy_records}**")
lines.append(f"- Docking shards: **{total_shards}**, failed shards: **{total_failed_shards}**")
lines.append(f"- Recorded failures across all stages: **{len(failure_rows)}**")
lines.append(f"- Case-insensitive input filename resolutions: **{len(case_match_rows)}**")
lines.append("")
lines.append("## 2. Search parameters")
lines.append("")
lines.append("| Parameter | Value |")
lines.append("|---|---:|")
for key in [
    "spacing_A", "padding_each_side_A", "center_method", "ga_run_total_per_pocket",
    "ga_num_evals_per_run", "ga_num_generations_per_run", "ga_population_size",
    "rmstol_A", "actual_shards_per_pocket", "planned_concurrent_docking_tasks",
]:
    if key in parameters:
        lines.append(f"| {key} | {parameters[key]} |")
lines.append("")
lines.append("## 3. Completion status")
lines.append("")
lines.append("| Status | Pocket count |")
lines.append("|---|---:|")
for status in sorted(status_counts):
    lines.append(f"| {status} | {status_counts[status]} |")
lines.append("")
lines.append("## 4. Per-pocket results")
lines.append("")
lines.append("| CAS | Receptor | Pocket | Status | Runs recovered/requested | Best | Median | Mean ± SD | Grid npts |")
lines.append("|---|---|---|---|---:|---:|---:|---:|---|")
for r in pocket_rows:
    mean_sd = "NA" if not r["mean_energy_kcal_mol"] else f"{r['mean_energy_kcal_mol']} ± {r['sd_energy_kcal_mol']}"
    best = r["best_energy_kcal_mol"] or "NA"
    median = r["median_energy_kcal_mol"] or "NA"
    lines.append(
        f"| {r['cas']} | {r['receptor']} | {r['pocket']} | {r['status']} | "
        f"{r['energy_records_found']}/{r['requested_lga_runs']} | {best} | {median} | {mean_sd} | {r['npts']} |"
    )
lines.append("")
lines.append("## 5. Overall ranking by lowest predicted binding energy")
lines.append("")
lines.append("| Rank | CAS | Receptor | Pocket | Best energy (kcal/mol) | Median | Completion | Best DLG |")
lines.append("|---:|---|---|---|---:|---:|---:|---|")
for rank, r in enumerate(ranked, start=1):
    lines.append(
        f"| {rank} | {r['cas']} | {r['receptor']} | {r['pocket']} | {r['best_energy_kcal_mol']} | "
        f"{r['median_energy_kcal_mol']} | {r['completion_percent']}% | `{r['best_dlg']}` |"
    )
lines.append("")
lines.append("## 6. Interpretation boundaries")
lines.append("")
lines.append("- AutoDock4 scores are model-dependent predicted interaction energies; they are not experimental binding free energies and do not demonstrate receptor activation.")
lines.append(f"- The {ga_run_total} LGA runs are divided into independent parallel shards to use the CPU allocation. Energy distributions are aggregated globally, but AutoDock's built-in RMSD clustering remains shard-local rather than one global {ga_run_total}-pose cluster analysis.")
lines.append("- Primary comparisons should therefore use completion-checked energy distributions, recurrent poses/interactions, and receptor-pocket consistency—not a single lowest score alone.")
lines.append("")
lines.append("## 7. Output files")
lines.append("")
for name, description in [
    ("pocket_summary.tsv", "one-row-per-pocket statistical summary"),
    ("top_hits.tsv", "all pockets ranked by lowest predicted energy"),
    ("all_binding_energies.tsv", "every recovered energy value with seeds and source DLG"),
    ("docking_shard_results.tsv", "execution status for every parallel shard"),
    ("failures.tsv", "all input, preparation, grid and docking failures"),
    ("job_manifest.tsv", "input paths, grid center, npts and parameter directories"),
    ("receptor_md5.tsv", "checksum consistency of complete receptor copies"),
    ("run_parameters.tsv", "complete run parameter record"),
    ("case_insensitive_matches.tsv", "audit trail of expected filenames resolved by case-insensitive matching"),
]:
    lines.append(f"- `{name}`: {description}.")

report_md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

plain = []
plain.append("AUTODOCK4 HIGH-PRECISION BATCH PIPELINE — MASTER REPORT")
plain.append("=" * 78)
plain.append(f"Parameter ID                : {parameter_id}")
plain.append(f"Pocket combinations         : {len(pocket_rows)}")
plain.append(f"Total requested LGA runs    : {total_requested}")
plain.append(f"Energy records recovered    : {total_energy_records}")
plain.append(f"Total docking shards        : {total_shards}")
plain.append(f"Failed docking shards       : {total_failed_shards}")
plain.append(f"Failures across all stages  : {len(failure_rows)}")
plain.append(f"Case-insensitive resolutions: {len(case_match_rows)}")
plain.append("")
plain.append("PER-POCKET RESULTS")
plain.append("-" * 78)
for r in pocket_rows:
    plain.append(
        f"{r['cas']} | {r['receptor']} | {r['pocket']} | {r['status']} | "
        f"runs {r['energy_records_found']}/{r['requested_lga_runs']} | "
        f"best {r['best_energy_kcal_mol'] or 'NA'} | median {r['median_energy_kcal_mol'] or 'NA'} | "
        f"mean {r['mean_energy_kcal_mol'] or 'NA'} | SD {r['sd_energy_kcal_mol'] or 'NA'}"
    )
plain.append("")
plain.append("RANKING")
plain.append("-" * 78)
for rank, r in enumerate(ranked, start=1):
    plain.append(
        f"{rank:>3}. {r['cas']} | {r['receptor']} | {r['pocket']} | "
        f"best={r['best_energy_kcal_mol']} kcal/mol | median={r['median_energy_kcal_mol']} | "
        f"completion={r['completion_percent']}%"
    )
plain.append("")
plain.append("CAUTION: docking scores are comparative model outputs, not proof of receptor activation or experimental affinity.")
plain.append("CAUTION: parallel sharding preserves the total number of LGA searches, but RMSD clustering is shard-local.")
report_txt_path.write_text("\n".join(plain) + "\n", encoding="utf-8")

readme_path.write_text(
    "Result priority:\n"
    "1. MASTER_REPORT.md / MASTER_REPORT.txt: human-readable overview.\n"
    "2. pocket_summary.tsv: principal statistical result table.\n"
    "3. all_binding_energies.tsv: full run-level evidence.\n"
    "4. top_hits.tsv: ranking only; do not interpret the minimum score alone.\n"
    "5. failures.tsv and docking_shard_results.tsv: completion and quality control.\n"
    "6. case_insensitive_matches.tsv: every filename accepted after a case-insensitive check.\n",
    encoding="utf-8",
)
PYEOF
report_rc=$?
if [ "$report_rc" -ne 0 ]; then
    log "ERROR: report aggregation failed with exit code $report_rc."
    log "Raw docking status table : $SHARD_RESULTS"
    log "Expected report path     : $RUN_DIR/MASTER_REPORT.txt"
    log "Expected Markdown report : $RUN_DIR/MASTER_REPORT.md"
    exit 5
fi

# Stable pointer to the newest successfully aggregated report run.
ln -sfn "$RUN_DIR" "$RESULT_ROOT/latest" 2>/dev/null || true

# =============================================================================
# 8. Final terminal summary
# =============================================================================
complete_pockets="$(awk -F '\t' 'NR>1 && $4=="COMPLETE" {n++} END{print n+0}' "$RUN_DIR/pocket_summary.tsv" 2>/dev/null || echo 0)"
partial_pockets="$(awk -F '\t' 'NR>1 && $4=="PARTIAL" {n++} END{print n+0}' "$RUN_DIR/pocket_summary.tsv" 2>/dev/null || echo 0)"
failed_pockets="$(awk -F '\t' 'NR>1 && $4=="FAILED" {n++} END{print n+0}' "$RUN_DIR/pocket_summary.tsv" 2>/dev/null || echo 0)"
failure_rows_count="$(awk 'END{print (NR>0 ? NR-1 : 0)}' "$FAILURES")"
case_match_count="$(awk 'END{print (NR>0 ? NR-1 : 0)}' "$CASE_MATCHES")"

log "============================================================"
log "Pipeline finished"
log "CAS directories processed : $cas_seen"
log "Receptor folders processed: $receptor_seen"
log "Pocket folders discovered : $pocket_seen"
log "Configurations generated  : $configured_count"
log "Complete pockets           : $complete_pockets"
log "Partial pockets            : $partial_pockets"
log "Failed pockets             : $failed_pockets"
log "Failure records            : $failure_rows_count"
log "Case-insensitive resolutions: $case_match_count"
log "Master report (text)       : $RUN_DIR/MASTER_REPORT.txt"
log "Master report (Markdown)   : $RUN_DIR/MASTER_REPORT.md"
log "Pocket summary             : $RUN_DIR/pocket_summary.tsv"
log "All energy values          : $RUN_DIR/all_binding_energies.tsv"
log "Overall ranking            : $RUN_DIR/top_hits.tsv"
log "Case-match audit           : $CASE_MATCHES"
log "Latest-run link            : $RESULT_ROOT/latest"
log "============================================================"

if [ "$RUN_DOCKING" = "1" ] && [ "$complete_pockets" -eq 0 ] && [ "$partial_pockets" -eq 0 ]; then
    log "ERROR: no docking result was recovered."
    exit 4
fi
exit 0
