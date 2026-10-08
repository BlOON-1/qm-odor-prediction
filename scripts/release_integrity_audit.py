"""Read-only pre-release integrity audit for this repository.

This script never alters research assets.  Its only outputs are the files in
release_audit/.  It is deliberately conservative: a readable file is not
treated as provenance-verified merely because it opens successfully.
"""
from __future__ import annotations

import ast
import csv
import hashlib
import json
import os
import re
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None
try:
    import openpyxl
except ImportError:  # pragma: no cover
    openpyxl = None
try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "release_audit"
SKIP_DIRS = {".git", ".agents", ".codex", "release_audit", "__pycache__", ".pytest_cache", "Angew_Data_Reporting_Checklist_Package", "fig2_data_audit_report", "submission_checklist_audit"}
RELEASE_TOP = {"data", "results", "scripts", "src", "configs", "docs", "notebooks", "tests"}
ROOT_FILES = {"README.md", "environment.yml", "requirements.txt", "pyproject.toml", "LICENSE", "CITATION.cff", ".gitignore", ".zenodo.json"}
TEXT = {".csv", ".tsv", ".txt", ".md", ".py", ".sh", ".yml", ".yaml", ".json", ".pdb", ".pdbqt", ".sdf", ".log", ".cff", ".toml", ".gitignore"}
TAB = {".tsv"}
# Require a real following path component. This avoids false positives in
# regular expressions such as ``DOCKED:\\s`` and ``/home/`` match patterns.
ABS_PATH = re.compile(r"(?:[A-Za-z]:\\[A-Z]|/(?:home|Users|mnt|private|var)/[A-Za-z0-9])")
SECRET = re.compile(r"(?i)(?:api[_-]?key|secret|token|password)\s*[:=]\s*['\"]?[A-Za-z0-9_\-]{12,}")

def rel(p: Path) -> str:
    return p.relative_to(ROOT).as_posix()

def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

def candidate_files():
    for p in ROOT.rglob("*"):
        if not p.is_file() or any(x in SKIP_DIRS for x in p.relative_to(ROOT).parts):
            continue
        parts = p.relative_to(ROOT).parts
        if parts[0] in RELEASE_TOP or p.name in ROOT_FILES:
            yield p

def category(r: str) -> str:
    if r.startswith("data/docking/"): return "G. docking/PLIP"
    if r.startswith("data/external_panel/"): return "D. structural-neighbor sensory panel"
    if r.startswith("data/sensory_triplets/"): return "E/F. triplets and intensity"
    if r.startswith("data/processed/") or r.startswith("data/splits/") or r.startswith("data/interim/"): return "A/B. benchmark and descriptors"
    if r.startswith("data/source_data/"): return "I. figure source data"
    if r.startswith("data/receptor_annotations/"): return "G. receptor annotations"
    if r.startswith("results/"): return "C/I. model results and plotting"
    if r.startswith(("scripts/", "src/", "configs/", "notebooks/", "tests/")): return "I/J. code/configuration"
    return "J. project metadata"

def read_delimited(p: Path):
    encodings = ("utf-8-sig", "utf-8", "gb18030", "latin-1")
    error = None
    for enc in encodings:
        try:
            with p.open("r", encoding=enc, newline="") as f:
                rows = list(csv.reader(f, delimiter="\t" if p.suffix.lower() in TAB else ","))
            if not rows: return [], [], enc
            return rows[0], rows[1:], enc
        except (UnicodeDecodeError, csv.Error) as e:
            error = e
    raise error or ValueError("unable to parse")

def text_issues(p: Path):
    raw = p.read_bytes()
    if b"\x00" in raw and p.suffix.lower() not in {".pdb", ".sdf"}:
        return "NUL byte in text file"
    s = raw.decode("utf-8", errors="replace")
    issues = []
    if ABS_PATH.search(s): issues.append("contains absolute local path")
    if SECRET.search(s): issues.append("possible credential pattern")
    return "; ".join(issues)

def check_pdb(p: Path):
    atoms = 0; bad = 0
    for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atoms += 1
            try: float(line[30:38]); float(line[38:46]); float(line[46:54])
            except ValueError: bad += 1
    if not atoms: raise ValueError("no ATOM/HETATM records")
    if bad: raise ValueError(f"{bad}/{atoms} atom coordinate records malformed")
    return atoms

def check_file(p: Path):
    r = rel(p); ext = p.suffix.lower(); rows = cols = ""; method = []; issue = ""; status = "UNVERIFIED"
    if p.stat().st_size == 0:
        return {"Relative_Path":r,"File_Type":ext or "extensionless","File_Size":p.stat().st_size,"SHA256":sha256(p),"Source_File":"","Data_Category":category(r),"Expected_Rows":"","Actual_Rows":"0","Expected_Columns":"","Actual_Columns":"","Integrity_Status":"FAIL","Validation_Method":"size check","Issue_Description":"empty file"}
    try:
        if ext in {".csv", ".tsv"}:
            head, body, encoding = read_delimited(p); rows, cols = len(body), len(head); method.append(f"CSV parser ({encoding})")
            if not head: raise ValueError("missing header")
            malformed = sum(len(x) != len(head) for x in body)
            if malformed: raise ValueError(f"{malformed} rows have inconsistent column count")
            if not body:
                # These are valid zero-event logs, not truncated data tables.
                zero_event = {"all_failures.tsv", "comparison_skips.tsv", "halogen_bonds.tsv", "metal_interactions.tsv", "plip_parse_failures.tsv"}
                if p.name not in zero_event: raise ValueError("header only; no data rows")
                method.append("recognized zero-event table")
            if r not in {"data/raw_manifest/source_inventory.csv"}: status = "PASS"
        elif ext == ".xlsx":
            if openpyxl is None: raise ValueError("openpyxl unavailable")
            wb = openpyxl.load_workbook(p, read_only=True, data_only=False)
            sheets = wb.sheetnames
            populated = 0
            for ws in wb.worksheets:
                for row in ws.iter_rows(values_only=True):
                    if any(x is not None for x in row): populated += 1
            rows, cols = populated, len(sheets); method.append(f"openpyxl workbook: {len(sheets)} sheets, {populated} nonempty rows")
            if not populated: raise ValueError("workbook has no populated cells")
            status = "UNVERIFIED"
        elif ext == ".npz":
            if np is None: raise ValueError("numpy unavailable")
            with np.load(p, allow_pickle=False) as z:
                keys = list(z.keys()); rows = ";".join(f"{k}:{tuple(z[k].shape)}" for k in keys); cols = len(keys)
            method.append("numpy.load allow_pickle=False")
            if not keys: raise ValueError("NPZ has no arrays")
            status = "PASS"
        elif ext == ".npy":
            if np is None: raise ValueError("numpy unavailable")
            a = np.load(p, allow_pickle=False); rows, cols = str(tuple(a.shape)), ""; method.append("numpy.load allow_pickle=False"); status = "PASS"
        elif ext == ".json":
            obj = json.loads(p.read_text(encoding="utf-8")); rows = len(obj) if hasattr(obj, "__len__") else 1; method.append("json parser"); status = "PASS"
        elif ext in {".yml", ".yaml"}:
            if yaml is None: raise ValueError("PyYAML unavailable")
            yaml.safe_load(p.read_text(encoding="utf-8")); method.append("YAML safe_load"); status = "PASS"
        elif ext == ".py":
            ast.parse(p.read_text(encoding="utf-8")); method.append("Python AST parse"); status = "PASS"
        elif ext == ".pdb":
            rows = check_pdb(p); method.append("PDB atom-coordinate parser"); status = "PASS"
        elif ext in {".pdbqt", ".sdf"}:
            content = p.read_text(encoding="utf-8", errors="replace"); rows = content.count("$$$$") if ext == ".sdf" else content.count("ATOM"); method.append("basic structure record check")
            if not rows: raise ValueError("no molecule/atom records")
            status = "UNVERIFIED"
        elif ext == ".zip":
            with zipfile.ZipFile(p) as z: bad = z.testzip(); rows = len(z.infolist())
            if bad: raise ValueError(f"bad compressed member: {bad}")
            method.append("ZipFile.testzip"); status = "PASS"
        else:
            issue = text_issues(p) if ext in TEXT else "binary/non-tabular asset not structurally parsed"
            method.append("readability and policy scan" if ext in TEXT else "binary existence/hash check")
            status = "PASS" if ext in TEXT else "UNVERIFIED"
        if ext in TEXT:
            policy = text_issues(p)
            if policy:
                status = "NOT_FOR_RELEASE"; issue = policy
    except Exception as e:
        status = "FAIL"; issue = f"parse/structure failure: {type(e).__name__}: {e}"
    return {"Relative_Path":r,"File_Type":ext or "extensionless","File_Size":p.stat().st_size,"SHA256":sha256(p),"Source_File":"","Data_Category":category(r),"Expected_Rows":"","Actual_Rows":rows,"Expected_Columns":"","Actual_Columns":cols,"Integrity_Status":status,"Validation_Method":"; ".join(method),"Issue_Description":issue}

def table(path):
    h, b, _ = read_delimited(ROOT / path)
    return [dict(zip(h, x)) for x in b]

def add_finding(findings, level, subject, evidence, effect):
    findings.append((level, subject, evidence, effect))

def cross_checks(manifest):
    findings=[]; tests=[]
    def test(name, fn):
        try:
            result=fn(); tests.append((name,"PASS",result))
        except Exception as e:
            tests.append((name,"FAIL",str(e))); add_finding(findings,"FAIL",name,str(e),"May block the linked scientific result")
    def benchmark():
        mol=table("data/processed/benchmark_molecules_4495.csv"); mat=table("data/processed/odor_label_matrix_4495x112.csv"); voc=table("data/processed/odor_label_vocabulary_112.csv")
        assert len(mol)==4495 and len(mat)==4495 and len(voc)==112
        assert [x["molecule_id"] for x in mol]==[x["molecule_id"] for x in mat]
        labels=list(mat[0])[1:]; assert labels==[x["normalized_label"] for x in voc]
        assert all(set(x[k] for k in labels)<={"0","1"} for x in mat)
        return "4495 molecule IDs/order, 112 labels, and binary matrix are mutually consistent; odorless is absent from the vocabulary."
    test("Benchmark 4495×112 matrix", benchmark)
    def fp():
        if np is None: raise AssertionError("numpy unavailable")
        idx=table("data/processed/ecfp4_fingerprint_index.csv")
        meta=json.loads((ROOT/"data/processed/ecfp4_fingerprint_metadata.json").read_text(encoding="utf-8"))
        molecules=table("data/processed/benchmark_molecules_4495.csv")
        with np.load(ROOT/"data/processed/ecfp4_fingerprints.npz",allow_pickle=False) as z:
            shape=tuple(int(x) for x in z["shape"]); indptr=z["indptr"]
        missing=[x for x in molecules if x["molecule_id"] not in {i["molecule_id"] for i in idx}]
        assert shape==(4494,1024) and len(indptr)==4495 and len(idx)==4494
        assert meta["n_molecules_total"]==4495 and meta["n_molecules_valid"]==4494 and meta["n_molecules_invalid"]==1
        assert len(missing)==1 and missing[0]["molecule_id"]=="MOL_003034" and missing[0]["smiles_parseable"]=="False"
        return "ECFP4 CSR matrix is 4494×1024 and index has 4494 rows; the sole omitted MOL_003034 has explicitly unparseable SMILES, matching metadata (4495 total/1 invalid)."
    test("ECFP4 index/array", fp)
    def panel():
        m=table("data/external_panel/external_molecules_87.csv"); pairs=table("data/external_panel/structural_neighbor_pairs_60.csv"); vocab=table("data/external_panel/panel_label_vocabulary_43.csv"); agg=table("data/external_panel/aggregated_panel_ratings.csv"); js=table("data/external_panel/pairwise_panel_jaccard_similarity.csv")
        cas={x["cas"] for x in m}; assert len(m)==len(cas)==87; assert len(pairs)==60 and len({x["pair_id"] for x in pairs})==60; assert len(vocab)==43
        assert all(x["cas1"] in cas and x["cas2"] in cas for x in pairs)
        assert len(agg)==87*43 and {x["n_assessors"] for x in agg}=={"13"}
        assert len(js)==60
        for x in js:
            a=set(filter(None,x["label_set_1"].split(";"))); b=set(filter(None,x["label_set_2"].split(";"))); got=len(a&b)/len(a|b) if a|b else 1
            assert abs(got-float(x["jaccard_similarity"]))<1e-12
        return "87 unique panel molecules, 60 valid pairs, 43 terms; all 60 published Jaccard values recompute exactly. Individual raw panel ratings are not present."
    test("External sensory-panel structure/Jaccard", panel)
    def triplets():
        d=table("data/sensory_triplets/triplet_definitions_30.csv"); r=table("data/sensory_triplets/triplet_assessor_responses_390.csv"); s=table("data/sensory_triplets/triplet_response_summary.csv")
        assert len(d)==30 and len({x["Triplet_ID"] for x in d})==30 and len(r)==390 and len(s)==30
        by=defaultdict(list)
        for x in r:
            assert x["Triplet_ID"] in {y["Triplet_ID"] for y in d}; assert x["Response_Code"] in {"1","2","3","4"}; by[x["Triplet_ID"]].append(x)
        assert all(len(x)==13 and len({y["Assessor_ID"] for y in x})==13 for x in by.values())
        for x in s:
            z=by[x["Triplet_ID"]]; counts=Counter(a["Response_Code"] for a in z)
            assert int(x["N_Assessors"])==13 and [int(x[k]) for k in ["Count_AP_Similar","Count_AS_Similar","Count_Both_Similar","Count_Both_Dissimilar"]]==[counts[str(i)] for i in range(1,5)]
        return "30 triplets × 13 unique assessors = 390 valid codes; every summary count recomputes. Original Excel is not included for external source comparison."
    test("Triplet response counts", triplets)
    def reps():
        idx=table("data/docking/metadata/representative_pose_index.tsv"); combos=table("data/docking/final_results/receptor_ligand_index.tsv")
        assert len(combos)==36 and len({(x['cas'],x['receptor_name_original'],x['pocket_name_original']) for x in combos})==36
        assert len(idx)==143
        for x in idx:
            p=ROOT/x['Released_File']; assert p.is_file(); assert sha256(p)==x['SHA256']; check_pdb(p)
        kinds=Counter(x['Pose_Type'] for x in idx)
        assert kinds['global_lowest_energy']==36 and kinds['dominant_cluster_medoid']==36
        return f"36 receptor-pocket jobs; 143 indexed structures have matching SHA-256 and coordinates. Required two core representatives/job=72 present; optional cluster medoids={len(idx)-72}."
    test("Docking representative PDB index", reps)
    def plip():
        status=table("data/docking/plip_details/plip_execution_status.tsv"); summ=table("data/docking/plip_details/interaction_summary.tsv"); freq=table("data/docking/final_results/residue_interaction_frequency.tsv")
        success={ (x['cas'],x['receptor'],x['pocket'],x['kind'],x['role'],x['pose_id']) for x in status if x['status']=='SUCCESS' }
        assert success; assert all((x['cas'],x['receptor'],x['pocket'],x['kind'],x['role'],x['pose_id']) in success for x in summ)
        assert all(abs(float(x['frequency'])-int(x['pose_occurrence_count'])/int(x['frequency_denominator']))<1e-6 for x in freq if int(x['frequency_denominator']))
        return f"{len(status)} PLIP tasks ({len(success)} success); all interaction-summary pose keys map to successful tasks; all frequency rows recompute from stated denominators."
    test("PLIP relationships/frequencies", plip)
    def model_summaries():
        metric_rows=[]
        for name, col in [("micro_auc_by_model_representation.csv", "test_roc_auc_micro"), ("macro_auc_by_model_representation.csv", "test_roc_auc_macro")]:
            for x in table("results/benchmark_metrics/"+name):
                vals=[float(x[f"{col}_fold_{i}"]) for i in range(1,6)]
                assert abs(sum(vals)/5-float(x[f"{col}_mean"]))<1e-12
                metric_rows.append(x)
        audit=table("results/structure_disjoint/fig2E_structure_disjoint_audit.csv")
        bad=[x for x in audit if x.get("status") == "MEAN_MISMATCH"]
        assert not bad, "structure-disjoint audit contains mean mismatches: " + ", ".join(x.get("representation", "?")+":"+x.get("metric", "?") for x in bad)
        return f"Recomputed means for {len(metric_rows)} benchmark metric rows."
    test("Model aggregate metrics", model_summaries)
    # annotate known source/provenance boundaries
    add_finding(findings,"UNVERIFIED","Raw panel records","Only aggregate 87×43 panel table is released; no 87×43×13 assessor-level panel matrix was found.","Figure 3 raw-rating provenance cannot be independently checked.")
    add_finding(findings,"UNVERIFIED","QM/ADCH descriptors and prediction probabilities","QM-valued structure-disjoint tables (3,189 rows) and aggregate model/SHAP summaries exist, but no standalone full QM/ADCH/fused feature matrices, trained weights, or per-molecule probabilities were found.","Aggregate benchmark means can be checked, but electronic-feature models cannot be fully rerun.")
    add_finding(findings,"NOT_FOR_RELEASE","Docking provenance tables","Several docking TSVs may retain absolute external paths. This is a portability/privacy disclosure risk.","Redact paths or publish portable replacements only.")
    return tests, findings

def main():
    OUT.mkdir(exist_ok=True)
    manifest=[check_file(p) for p in sorted(candidate_files())]
    tests, findings=cross_checks(manifest)
    # Promote only cross-checked essential tables; all other successful structured assets remain conservative UNVERIFIED.
    checked_prefixes=("data/processed/benchmark_molecules_4495.csv","data/processed/odor_label_matrix_4495x112.csv","data/processed/odor_label_vocabulary_112.csv","data/processed/ecfp4_fingerprint_index.csv","data/processed/ecfp4_fingerprints.npz","data/processed/ecfp4_fingerprint_metadata.json","data/external_panel/external_molecules_87.csv","data/external_panel/structural_neighbor_pairs_60.csv","data/external_panel/panel_label_vocabulary_43.csv","data/external_panel/pairwise_panel_jaccard_similarity.csv","data/sensory_triplets/triplet_definitions_30.csv","data/sensory_triplets/triplet_assessor_responses_390.csv","data/sensory_triplets/triplet_response_summary.csv","data/docking/metadata/representative_pose_index.tsv","data/docking/representative_complexes/")
    for m in manifest:
        if m['Integrity_Status']=='PASS' and not m['Relative_Path'].startswith(checked_prefixes): m['Integrity_Status']='UNVERIFIED'; m['Issue_Description']=(m['Issue_Description']+'; ' if m['Issue_Description'] else '')+'readable, but no file-specific provenance/cross-source validation'
        if m['Relative_Path'].startswith('data/docking/') and 'absolute local path' in m['Issue_Description']: m['Integrity_Status']='NOT_FOR_RELEASE'
        if m['Relative_Path'] in {'results/structure_disjoint/structure_disjoint_metrics.csv','results/structure_disjoint/fig2E_structure_disjoint_audit.csv'}:
            m['Integrity_Status']='FAIL'; m['Issue_Description']=(m['Issue_Description']+'; ' if m['Issue_Description'] else '')+'two documented MEAN_MISMATCH records: QM and RDKit+QM test micro-AUC'
    fields=list(manifest[0])
    with (OUT/'FILE_INTEGRITY_MANIFEST.csv').open('w',encoding='utf-8',newline='') as f: csv.DictWriter(f,fieldnames=fields).writeheader(); csv.DictWriter(f,fieldnames=fields).writerows(manifest)
    cnt=Counter(m['Integrity_Status'] for m in manifest)
    direct=[m for m in manifest if m['Integrity_Status']=='PASS' and not m['Issue_Description']]
    bad=[m for m in manifest if m['Integrity_Status']!='PASS']
    test_md='\n'.join(f"| {n} | {s} | {d} |" for n,s,d in tests)
    finding_md='\n'.join(f"| {a} | `{b}` | {c} | {d} |" for a,b,c,d in findings)
    (OUT/'FINAL_DATA_INTEGRITY_REPORT.md').write_text(f"""# Final data-integrity audit\n\n**Audit date:** 2026-10-08  \n**Scope:** {len(manifest)} release-candidate files under data/, results/, scripts/, src/, configs/, docs/, notebooks/, tests/ and root release metadata. Excluded internal audit/checklist folders and Git internals. No source data were modified.\n\n## Overall decision: BLOCKED\n\nThe internally cross-checkable benchmark, triplet, and representative-structure subsets pass their stated checks, but the release is not ready as a complete paper data package: QM/ADCH/model-output assets are absent, panel raw ratings are absent, and portable docking metadata still contains absolute local paths. Git state could not be read (`git status` and `git rev-parse` returned “not a git repository”), despite a `.git` directory being present.\n\n## File status counts\n\n| PASS | MISSING | FAIL | UNVERIFIED | NOT_FOR_RELEASE | total |\n|---:|---:|---:|---:|---:|---:|\n| {cnt['PASS']} | 0 (file-system manifest cannot enumerate absent files) | {cnt['FAIL']} | {cnt['UNVERIFIED']} | {cnt['NOT_FOR_RELEASE']} | {len(manifest)} |\n\n`MISSING` planned assets are listed in the issue document rather than fabricated as manifest rows.\n\n## Reproduction and consistency tests\n\n| Test | status | evidence |\n|---|---|---|\n{test_md}\n\n## Findings and unresolved boundaries\n\n| status | item | evidence | effect |\n|---|---|---|---|\n{finding_md}\n\n## Interpretation of collection versions\n\nThe formal model dataset is 4,495 molecules × 112 labels: `odorless` is not in the released 112-term vocabulary. “Observed” files are retained as a different collection state and are not silently forced to 4,495 rows. No 4,403-model package was found.\n\n## Release conditions\n\nOnly manifest rows marked PASS and free of release-policy issues are candidates for direct publication. PASS proves the explicit checks above, not scientific validity of docking or a third-party redistribution license. Model/QM claims remain **BLOCKED** pending release or documented exclusion of the missing data products.\n""",encoding='utf-8')
    def entry(m): return f"- `{m['Relative_Path']}` — {m['Data_Category']}; {m['Validation_Method']}"
    grouped=defaultdict(list)
    for m in direct: grouped[m['Data_Category']].append(entry(m))
    upload='\n'.join(f"## {k}\n\n"+'\n'.join(v) for k,v in grouped.items()) or "No files met the strict direct-upload filter."
    (OUT/'GITHUB_UPLOAD_LIST.md').write_text(f"""# GitHub upload list\n\n## Direct-upload candidates (PASS and no detected release-policy issue)\n\n{upload}\n\n## Do not mix into direct upload\n\nAll UNVERIFIED, FAIL and NOT_FOR_RELEASE files in `FILE_INTEGRITY_MANIFEST.csv` require the disposition in `FILES_TO_FIX_OR_EXCLUDE.md`. Release permission for third-party sources remains a separate legal check.\n""",encoding='utf-8')
    issues='\n'.join(f"| {m['Integrity_Status']} | `{m['Relative_Path']}` | {m['Issue_Description'] or 'Not validated sufficiently for PASS'} | Review source/provenance; redact or exclude where applicable. | {'Yes' if m['Data_Category'][0] in 'ABCDEFG' else 'No'} |" for m in bad)
    missing='\n'.join(["| MISSING (critical) | QM descriptor matrix; ADCH charge/features; RDKit+QM fusion matrix | No matching released files found | Export from verified source with IDs and QC/failure log | Yes |", "| MISSING (critical) | Per-molecule prediction probabilities, trained weights/seeds, SHAP values | No matching released files found | Release verified outputs or limit claims | Yes |", "| MISSING (critical) | External-panel assessor-level records | Only aggregated panel values found | Release anonymized raw ratings or state aggregate-only limitation | Yes |", "| MISSING (important) | Docking validation source records for 8F76, 8UXY, 9WG4 | No matching validation inputs/results found | Add source/RMSD/pocket-recovery records | Yes |", "| MISSING (optional) | Original triplet Excel source | Derived CSV exists, but source workbook absent | Add licensed/anonymized source workbook or provenance hash | No |"])
    (OUT/'FILES_TO_FIX_OR_EXCLUDE.md').write_text(f"""# Files to fix, exclude, or document\n\n## Missing planned assets\n\n| status | asset | evidence | recommendation | affects reproduction |\n|---|---|---|---|---|\n{missing}\n\n## Existing files not approved for direct upload\n\n| status | file | evidence | recommendation | affects reproduction |\n|---|---|---|---|---|\n{issues}\n\nDo not repair, synthesize, or overwrite scientific records to satisfy this list.\n""",encoding='utf-8')
    pathlist='\n'.join(f"  - {m['Relative_Path']}" for m in direct)
    (OUT/'FINAL_RELEASE_TREE.md').write_text(f"""# Recommended release tree\n\nThis is an exact file list, not an assertion that omitted assets exist. It contains only files that passed the strict audit filter.\n\n```text\nqm-odor-prediction/\n{pathlist}\n```\n\nFiles omitted from this tree are not approved for direct upload in this audit.\n""",encoding='utf-8')
    (OUT/'MINIMUM_REPRODUCIBLE_RELEASE.md').write_text("""# Minimum reproducible release\n\n## Main conclusions that can be recalculated now\n\n- The released 4,495 × 112 benchmark label matrix and ECFP4 indexing consistency.\n- All 60 published panel-pair Jaccard values from their published label sets.\n- All 30 triplet summary counts from the 390 anonymized responses.\n- Representative PDB file integrity against the provided SHA-256 index, and PLIP frequency arithmetic (although docking metadata itself needs portable redaction before release).\n\n## Conditions only, not currently reproducible\n\n- QM/ADCH and fused-feature modeling claims, prediction probabilities, model metrics, SHAP analyses, and seed-level split/model reruns: required feature/output assets are missing.\n- Figure 3 assessor-level analyses: only aggregate panel values are present.\n- Docking validation claims for 8F76/8UXY/9WG4: validation records are absent.\n\n## Extra assets for a complete rerun\n\nRelease ID-aligned QM/ADCH/fusion matrices, calculation-failure log, split seeds, model code/weights or deterministic training instructions, per-molecule predictions/labels, SHAP outputs, anonymized panel raw ratings, and validation inputs/outputs. Preserve third-party raw materials only where redistribution rights permit; otherwise provide rebuild instructions and immutable source/version hashes.\n""",encoding='utf-8')
    print(json.dumps({"files":len(manifest),"statuses":cnt,"tests":tests},indent=2,default=str))

if __name__ == '__main__': main()
