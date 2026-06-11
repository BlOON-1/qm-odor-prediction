from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_cas_set(csv_path: Path) -> set[str]:
    with csv_path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or "cas" not in reader.fieldnames:
            raise ValueError(f"在文件中找不到 'cas' 列: {csv_path}")
        cas_set: set[str] = set()
        for row in reader:
            cas = (row.get("cas") or "").strip()
            if cas:
                cas_set.add(cas)
    return cas_set


def read_rows(csv_path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with csv_path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"空表头: {csv_path}")
        rows = list(reader)
    return list(reader.fieldnames), rows


def write_rows(csv_path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="按各 fold 的 cas 列表切分 ECFP4_data.csv，并输出到 ECFP4/ 目录。"
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="数据根目录（默认：脚本所在目录）",
    )
    parser.add_argument(
        "--ecfp4-file",
        type=Path,
        default=None,
        help="ECFP4 数据文件（默认：<base-dir>/ECFP4_data.csv）",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="输出目录（默认：<base-dir>/ECFP4）",
    )
    args = parser.parse_args()

    base_dir: Path = args.base_dir.resolve()
    ecfp4_file: Path = (args.ecfp4_file or (base_dir / "ECFP4_data.csv")).resolve()
    out_dir: Path = (args.out_dir or (base_dir / "ECFP4")).resolve()

    ecfp4_fields, ecfp4_rows = read_rows(ecfp4_file)
    if "cas" not in ecfp4_fields:
        raise ValueError(f"ECFP4 文件中找不到 'cas' 列: {ecfp4_file}")

    data_cas_all = {(r.get("cas") or "").strip() for r in ecfp4_rows if (r.get("cas") or "").strip()}

    fold_dirs = sorted([p for p in base_dir.glob("fold_*") if p.is_dir()], key=lambda p: p.name)
    if not fold_dirs:
        raise ValueError(f"未找到 fold_* 目录: {base_dir}")

    print(f"base_dir: {base_dir}")
    print(f"ecfp4_file: {ecfp4_file}")
    print(f"out_dir: {out_dir}")
    print(f"ECFP4 总行数: {len(ecfp4_rows)} | 唯一 cas 数: {len(data_cas_all)}")

    for fold_dir in fold_dirs:
        fold_name = fold_dir.name  # e.g. fold_01
        suffix = fold_name.split("_", 1)[-1]
        train_src = fold_dir / f"{fold_name}.csv"
        test_src = fold_dir / f"test_{suffix}.csv"

        if not train_src.exists():
            raise FileNotFoundError(f"缺少训练集文件: {train_src}")
        if not test_src.exists():
            raise FileNotFoundError(f"缺少测试集文件: {test_src}")

        train_cas = read_cas_set(train_src)
        test_cas = read_cas_set(test_src)

        overlap = train_cas & test_cas
        if overlap:
            # 以测试集优先：从训练集中移除重叠 cas
            train_cas = train_cas - overlap

        union_cas = train_cas | test_cas
        missing_in_ecfp4 = union_cas - data_cas_all
        unused_in_fold = data_cas_all - union_cas

        train_rows_out = [r for r in ecfp4_rows if (r.get("cas") or "").strip() in train_cas]
        test_rows_out = [r for r in ecfp4_rows if (r.get("cas") or "").strip() in test_cas]

        out_fold_dir = out_dir / fold_name
        train_out = out_fold_dir / f"{fold_name}.csv"
        test_out = out_fold_dir / f"test_{suffix}.csv"
        write_rows(train_out, ecfp4_fields, train_rows_out)
        write_rows(test_out, ecfp4_fields, test_rows_out)

        print("\n" + "-" * 60)
        print(f"{fold_name}")
        print(f"train cas: {len(train_cas)} | test cas: {len(test_cas)} | overlap cas: {len(overlap)}")
        print(f"输出行数 train: {len(train_rows_out)} -> {train_out}")
        print(f"输出行数 test : {len(test_rows_out)} -> {test_out}")
        if missing_in_ecfp4:
            print(f"警告：fold 中 {len(missing_in_ecfp4)} 个 cas 在 ECFP4_data.csv 中找不到（示例前 5 个）：{sorted(missing_in_ecfp4)[:5]}")
        if unused_in_fold:
            print(f"提示：ECFP4_data.csv 中 {len(unused_in_fold)} 个 cas 未被该 fold 的 train/test 覆盖（示例前 5 个）：{sorted(unused_in_fold)[:5]}")
        if overlap:
            print(f"警告：train/test 出现 {len(overlap)} 个重复 cas（已按“测试集优先”处理）（示例前 5 个）：{sorted(overlap)[:5]}")

    print("\n完成。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

