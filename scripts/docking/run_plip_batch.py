#!/usr/bin/env python3
"""Run PLIP XML analysis with bounded parallelism and content-aware resume."""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from xml.etree import ElementTree

from pipeline_common import analysis_run_dir, read_tsv, update_stage_status, write_tsv


def valid_xml(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size < 50: return False
    try: ElementTree.parse(str(path)); return True
    except Exception: return False


def execute(row, raw_root: Path, skip: bool):
    key = "{}__{}__{}__{}__{}".format(row["cas"], row["receptor"], row["pocket"], row["kind"], row["role"])
    output = raw_root / key
    output.mkdir(parents=True, exist_ok=True)
    xmls = sorted(output.glob("*.xml"))
    if skip and len(xmls) == 1 and valid_xml(xmls[0]):
        return dict(row, plip_xml=str(xmls[0]), status="SKIPPED_COMPLETE", exit_code=0, message="validated existing PLIP XML")
    for path in xmls:
        path.unlink()
    command = ["plip", "-f", row["complex_pdb"], "-x", "-o", str(output)]
    try:
        completed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except Exception as exc:
        return dict(row, plip_xml="", status="FAILED", exit_code=127, message="PLIP launch failed: {!r}".format(exc))
    (output / "plip_command.log").write_text(completed.stdout or "", encoding="utf-8", errors="replace")
    xmls = sorted(output.glob("*.xml"))
    if completed.returncode == 0 and len(xmls) == 1 and valid_xml(xmls[0]):
        status, message, xml = "SUCCESS", "PLIP XML generated", str(xmls[0])
    else:
        status, message, xml = "FAILED", "PLIP exit {}; {} valid XML files".format(completed.returncode, sum(valid_xml(p) for p in xmls)), ""
    return dict(row, plip_xml=xml, status=status, exit_code=completed.returncode, message=message)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-run", type=Path, default=None)
    parser.add_argument("--parallel", type=int, default=int(os.environ.get("PLIP_PARALLEL", "16")))
    args = parser.parse_args()
    analysis = args.analysis_run or analysis_run_dir()
    try:
        allocated_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", os.environ.get("POST_CPUS", "32")))
    except ValueError as exc:
        raise RuntimeError("SLURM_CPUS_PER_TASK/POST_CPUS must be an integer") from exc
    if allocated_cpus < 1:
        raise RuntimeError("allocated CPU count must be positive")
    workers = min(max(1, args.parallel), 16, allocated_cpus)
    complex_rows = read_tsv(analysis / "complexes" / "complex_inventory.tsv")
    rows = [r for r in complex_rows if r["status"] in {"SUCCESS", "SKIPPED_COMPLETE"}]
    raw_root = analysis / "interactions" / "raw_plip"
    raw_root.mkdir(parents=True, exist_ok=True)
    skip = os.environ.get("SKIP_EXISTING", "1") == "1"
    results = [dict(r, plip_xml="", status="FAILED", exit_code=125,
                    message="complex unavailable: " + r.get("message", ""))
               for r in complex_rows if r["status"] not in {"SUCCESS", "SKIPPED_COMPLETE"}]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(execute, row, raw_root, skip) for row in rows]
        for future in as_completed(futures): results.append(future.result())
    results.sort(key=lambda r: (r["cas"], r["receptor"], r["pocket"], r["kind"], r["role"]))
    fields = ["cas", "receptor", "pocket", "kind", "role", "pose_id", "complex_pdb", "plip_xml", "status", "exit_code", "message"]
    write_tsv(analysis / "interactions" / "plip_execution_status.tsv", results, fields)
    rep_fail = sum(r["kind"] == "representative" and r["status"] == "FAILED" for r in results)
    failures = sum(r["status"] == "FAILED" for r in results)
    status = "SUCCESS" if results and not failures else ("WARNING" if results and not rep_fail else "FAIL")
    update_stage_status(analysis, "plip", status, "{} successful/reused; {} failures; workers={}".format(len(results) - failures, failures, workers), 0 if results and not rep_fail else 1)
    return 0 if results and not rep_fail else 1


if __name__ == "__main__":
    raise SystemExit(main())
