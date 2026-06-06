"""
Module 4: 训练集扩充收集脚本。

从 GraphPPIS（Train_335 + Test_60 + UBtest_31，已剔除与本项目 Dset_186-local /
PDBtest_315-local 评测集及原有 500 复合物语料的重叠）筛出的候选 PDB+链 列表中，
下载生物组装体并用与主语料完全一致的 ΔSASA 管线（阈值 2.0 A^2）生成弱标签。

输出（与 benchmark 流程同构，可直接喂给 extract_esm_embeddings / build_multimodal_dataset）:
    data/processed/train_expansion_labels.csv
    data/processed/train_expansion_manifest.csv

用法:
    $env:PYTHONPATH="src"
    python scripts/collect_train_expansion.py --max-usable 200
"""

import argparse
import csv
import sys
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(errors="replace")

SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))
sys.path.insert(0, str(SCRIPTS_DIR.parent / "src"))

import run_benchmark_eval as rbe
from sasa_project.paths import PROCESSED_DATA_DIR, PDB_COMPLEX_DIR


def load_candidates(path):
    entries = []
    with open(path) as f:
        for row in csv.DictReader(f):
            entries.append((row["pdb_id"].strip().upper(), row["target_chain"].strip().upper()))
    return entries


def main():
    parser = argparse.ArgumentParser(description="收集训练集扩充数据（GraphPPIS 经典 PPI 复合物）。")
    parser.add_argument(
        "--candidates",
        default="data/raw/graphppis/train_expansion_ids.csv",
        help="候选 PDB+target_chain 列表 CSV（pdb_id,target_chain）。",
    )
    parser.add_argument(
        "--output-dir",
        default=str(PDB_COMPLEX_DIR / "train_expansion"),
        help="扩充结构文件存放目录。",
    )
    parser.add_argument("--threshold", type=float, default=2.0, help="Delta-SASA 界面阈值（与主语料一致）。")
    parser.add_argument("--max-usable", type=int, default=200, help="达到该可用复合物数量后停止。")
    args = parser.parse_args()

    candidates = load_candidates(args.candidates)
    print(f"候选条目: {len(candidates)}（目标可用: {args.max_usable}）")

    # 构造 benchmark 风格 entry（PDBID+CHAIN），复用下载逻辑
    entries = [f"{pdb}{chain}" for pdb, chain in candidates]
    items = rbe.download_benchmark_structures(entries, args.output_dir)
    print(f"成功获取结构: {len(items)} / {len(entries)}")

    all_rows = []
    manifest_rows = []
    usable = 0
    skipped_chain = 0
    for item in items:
        if usable >= args.max_usable:
            break
        rows = rbe.process_benchmark_complex(
            item["pdb_path"],
            delta_sasa_threshold=args.threshold,
            target_chain=item["target_chain"],
        )
        if not rows:
            skipped_chain += 1
            continue
        all_rows.extend(rows)
        first = rows[0]
        manifest_rows.append(
            {
                "entry_id": item["entry_id"],
                "pdb_id": first["pdb_id"],
                "target_chain": first["target_chain"],
                "partner_chain": first["partner_chain"],
                "pdb_path": str(item["pdb_path"]),
            }
        )
        usable += 1

    label_out = PROCESSED_DATA_DIR / "train_expansion_labels.csv"
    manifest_out = PROCESSED_DATA_DIR / "train_expansion_manifest.csv"
    rbe.save_benchmark_labels(all_rows, label_out)
    rbe.save_benchmark_manifest(manifest_rows, manifest_out)

    pos = sum(int(r["label"]) for r in all_rows)
    total = len(all_rows)
    print("\n=== 扩充集统计 ===")
    print(f"可用复合物: {usable}（因链数/链长过滤跳过: {skipped_chain}）")
    print(f"目标链残基: {total}")
    if total:
        print(f"正例(界面): {pos} ({pos/total*100:.1f}%)  负例: {total-pos} ({(total-pos)/total*100:.1f}%)")


if __name__ == "__main__":
    main()
