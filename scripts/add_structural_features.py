"""
快速为现有多模态数据集添加结构特征（二面角 + HSE + 疏水性）。

无需重新运行 ESM 提取，直接在现有 CSV 基础上追加列。

使用方法:
    python scripts/add_structural_features.py \
        --input  data/processed/multimodal_residue_dataset_650m.csv \
        --manifest data/processed/complex_manifest.csv \
        --output data/processed/multimodal_residue_dataset_650m_struct.csv
"""

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from sasa_project.residue_features import (
    extract_structural_features,
    residue_key,
    resolve_pdb_path,
)
from sasa_project.sasa import parse_pdb


STRUCT_COLS = ["sin_phi", "cos_phi", "sin_psi", "cos_psi", "hse_up", "hse_dn", "hydrophobicity"]
STRUCT_DEFAULTS = {
    "sin_phi": "0.00000000",
    "cos_phi": "1.00000000",
    "sin_psi": "0.00000000",
    "cos_psi": "1.00000000",
    "hse_up": "0.0",
    "hse_dn": "0.0",
    "hydrophobicity": "0.0000",
}


def load_manifest(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {
        (row["pdb_id"], row["target_chain"], row["partner_chain"]): row
        for row in rows
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    manifest = load_manifest(args.manifest)

    # 缓存：每个复合体只解析一次 PDB
    struct_cache = {}

    def get_struct_features(pdb_id, target_chain, partner_chain, chain_id, res_id, ins_code):
        complex_key = (pdb_id, target_chain, partner_chain)
        if complex_key not in struct_cache:
            manifest_row = manifest.get(complex_key)
            if manifest_row is None:
                struct_cache[complex_key] = {}
            else:
                try:
                    pdb_path = resolve_pdb_path(manifest_row["pdb_path"])
                    atoms = parse_pdb(str(pdb_path))
                    struct_cache[complex_key] = extract_structural_features(atoms, chain_id)
                except Exception as exc:
                    print(f"[warn] {pdb_id}: {exc}")
                    struct_cache[complex_key] = {}

        rkey = residue_key(chain_id, res_id, ins_code)
        return struct_cache[complex_key].get(rkey, {})

    with open(args.input, newline="") as fin:
        reader = csv.DictReader(fin)
        existing_cols = reader.fieldnames or []

        # 如果已有结构特征列，跳过已有的
        new_cols = [c for c in STRUCT_COLS if c not in existing_cols]
        if not new_cols:
            print("数据集已包含所有结构特征列，无需添加。")
            return

        # 确定插入位置：在 label_threshold 之后，esm_ 列之前
        insert_after = "label_threshold"
        if insert_after in existing_cols:
            insert_pos = existing_cols.index(insert_after) + 1
        else:
            insert_pos = len(existing_cols)

        out_cols = existing_cols[:insert_pos] + new_cols + existing_cols[insert_pos:]

        total = 0
        with open(args.output, "w", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=out_cols)
            writer.writeheader()

            for i, row in enumerate(reader):
                feats = get_struct_features(
                    row["pdb_id"],
                    row["target_chain"],
                    row["partner_chain"],
                    row["chain_id"],
                    row["residue_id"],
                    row.get("insertion_code", ""),
                )
                for col in new_cols:
                    val = feats.get(col)
                    if val is None:
                        row[col] = STRUCT_DEFAULTS[col]
                    elif col in ("hse_up", "hse_dn"):
                        row[col] = f"{float(val):.1f}"
                    elif col == "hydrophobicity":
                        row[col] = f"{float(val):.4f}"
                    else:
                        row[col] = f"{float(val):.8f}"
                writer.writerow(row)
                total += 1
                if (i + 1) % 10000 == 0:
                    print(f"  处理进度: {i+1} 行...")

    print(f"完成：{total} 行写入 {args.output}")
    print(f"新增列：{new_cols}")


if __name__ == "__main__":
    main()
