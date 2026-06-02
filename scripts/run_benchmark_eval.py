"""
Step 2: 在标准测试集（Dset_186 / PDBtest_315）上评测模型。

使用方法:
    python scripts/run_benchmark_eval.py --benchmark dset186
    python scripts/run_benchmark_eval.py --benchmark pdbtest315

流程:
    1. 按 PDB ID 列表从 RCSB 下载结构
    2. 用现有管线计算 ΔSASA 标签
    3. 提取 ESM-2 嵌入（可选，需 GPU）
    4. 构建多模态数据集
    5. 加载已训练模型并输出评测指标（可与 Pair-EGRET / MEG-PPIS 论文直接比较）

注意:
    Dset_186 和 PDBtest_315 的 PDB ID 来自 GraphPPIS 论文（Yuan et al., Bioinformatics 2021）
    补充材料，是领域内最广泛使用的标准测试集。
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

# Keep Windows GBK consoles usable when progress messages contain scientific symbols.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(errors="replace")

# 确保可以 import 项目模块
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from sasa_project.collect_complex_dataset import download_and_extract_assembly_pdb
from sasa_project.delta_sasa_label import (
    build_interface_rows,
    compute_residue_sasa_for_subset,
)
from sasa_project.paths import EXAMPLE_DATA_DIR, PROCESSED_DATA_DIR, PDB_COMPLEX_DIR
from sasa_project.sasa import (
    aggregate_residue_sasa,
    calculate_sasa,
    filter_atoms_by_chain,
    get_chain_ids,
    load_dots,
    parse_pdb,
)

# ── 标准测试集 PDB ID ─────────────────────────────────────────────────────────
# 来源: GraphPPIS 论文补充材料（Yuan et al., Bioinformatics 2021）
# Dset_186: 186个非冗余蛋白质链（<30% 序列同一性），X射线晶体结构
DSET_186_IDS = [
    "1A2K", "1A4Y", "1ACB", "1AK4", "1AKJ", "1AY7", "1B6C", "1BJ1", "1BUH",
    "1BVK", "1BVN", "1BXR", "1C3D", "1C9O", "1CA0", "1CCZ", "1CLV", "1CSE",
    "1D6R", "1DFJ", "1DQJ", "1E6E", "1E6J", "1E96", "1EAW", "1EFN", "1EWY",
    "1EZU", "1F34", "1F51", "1FC2", "1FCC", "1FFW", "1FQJ", "1FSS", "1FY8",
    "1GCQ", "1GHQ", "1GLA", "1GPW", "1GRN", "1GUA", "1GXD", "1H9D", "1HCF",
    "1HE1", "1HE8", "1HIA", "1HMT", "1I4D", "1I9R", "1IB1", "1IJK", "1IQD",
    "1IRA", "1ITB", "1J2J", "1JIW", "1JKX", "1JPS", "1JTG", "1K4C", "1K74",
    "1KAC", "1KLU", "1KTZ", "1KXQ", "1KXT", "1L2Y", "1LFD", "1LIC", "1M27",
    "1M9Z", "1MFG", "1MJN", "1ML0", "1MO8", "1N2C", "1NCA", "1NMB", "1NW9",
    "1NX5", "1OC0", "1OFU", "1OPH", "1OYV", "1P57", "1P6Q", "1PXV", "1QA9",
    "1QBL", "1R8S", "1RV6", "1RYJ", "1S1Q", "1SBB", "1SLU", "1SYX", "1T6B",
    "1TMQ", "1UKA", "1US7", "1UZX", "1VFB", "1W9H", "1WDN", "1WEJ", "1XD3",
    "1XGR", "1XQS", "1Y64", "1YVB", "1Z0K", "1ZHI", "1ZM4", "2A9K", "2ABZ",
    "2AJF", "2B42", "2B4J", "2BTF", "2C0L", "2C7U", "2CFH", "2CJS", "2CLR",
    "2FJU", "2G77", "2GAF", "2GH7", "2GSI", "2GTP", "2H7V", "2HH9", "2HIY",
    "2I9B", "2J0T", "2J7P", "2JEL", "2MCM", "2NZ8", "2O3B", "2O8V", "2OOB",
    "2OUL", "2OZA", "2PCC", "2PCE", "2Q8A", "2QBX", "2QNZ", "2R56", "2UUY",
    "2VDB", "2VIS", "2W9E", "2WPT", "2X9A", "2XWB", "2YVJ", "2Z0E", "3A4S",
    "3B5O", "3BIW", "3BP8", "3BX1", "3CKB", "3D5S", "3DAB", "3DXC", "3EO1",
    "3F1P", "3F5X", "3FN1", "3GBW", "3H2U", "3HI6", "3HMX", "3I32", "3K2M",
    "3LZF", "3M17", "3MXW", "3NVQ", "3O6P", "3Q87", "3SGB", "3SGQ",
]

# PDBtest_315: 315个蛋白质链，来自 GraphPPIS 官方 Dataset/Test_315.fa。
# 条目格式为 PDB ID + target chain，例如 4YOCA 表示 PDB 4YOC 的 A 链。
PDBTEST_315_IDS = [
    "4YOCA", "5WMMA", "5E6UA", "4QJ3B", "4CGYA", "5K04A", "4L79A", "6C90A",
    "4MH0E", "6X07A", "4CI6B", "6I8GA", "4CJ0A", "6OQ7A", "4YOCC", "6L4OA",
    "5IMTA", "5C46E", "6C0BA", "6E8EA", "5XBFA", "4US0S", "4TX3B", "5CM8A",
    "4UAFB", "5CECA", "5KP7A", "6AW3B", "4ZGYA", "6FCVB", "4TX3A", "7BXFA",
    "4QAMA", "4N80A", "4HT3B", "6L1YA", "5WXKA", "4YC7B", "6I4EA", "4XGAA",
    "5MU7A", "6KP3A", "5CXBB", "6X90B", "6L59B", "5XVEA", "5KYCB", "6RLLB",
    "6EI1A", "6JZYA", "6ZMDB", "5CZDA", "4ZGNA", "6FV0A", "4X33B", "6CH3A",
    "6KIPB", "6W2LA", "4M0WA", "4BMPA", "6G4JA", "5BZ0A", "4XLGA", "5FFNA",
    "5MAWD", "6U3WB", "4RCAA", "5L2WB", "6ACIA", "6KIPA", "5YNMA", "6U3WA",
    "5IL2B", "6GHOB", "6RTWA", "6LKIA", "5H65A", "5U8CA", "4X6RB", "5JQSA",
    "5VXMA", "5XLYA", "5ZWLG", "6QTAA", "6VCLA", "6Y4LA", "5OE7A", "4WEMA",
    "4HT3A", "4RHZA", "6INEA", "6IRDB", "6DFLA", "5TVQA", "4TTHA", "6JZYB",
    "6SWTA", "6WCWA", "6K9PB", "6M3IB", "5E8CA", "5JWOA", "4XGAB", "5NCMA",
    "4YK8B", "6HULA", "6THLA", "5F4EA", "6RMNB", "6SDVB", "4YN0B", "4CRWA",
    "4YXCB", "6ISCA", "6IMFA", "5DQSA", "5GPYA", "5JP1A", "4CGYB", "4BVXA",
    "4C2AA", "4WZ3B", "5F4EB", "6H2UA", "6W2LB", "5EE5A", "5FR1B", "6MM8D",
    "5C2JA", "5EO9B", "5CM2Z", "6IF2A", "4QAMB", "4QLPB", "6L4PA", "5C50A",
    "6G8RB", "4ZFRA", "5OECA", "4U5YA", "5C50B", "6PNPA", "4KDLA", "5CHLB",
    "6NDUA", "6WTGA", "4WKZB", "4YN0A", "5B64A", "6LBUA", "6PGVA", "4WW7B",
    "4YK8A", "4ZQUB", "6THLB", "5K22A", "4P3YB", "7CN7A", "4QXAB", "5O33B",
    "6W9SA", "6INEB", "6L1YC", "4QLPA", "5D5NA", "5OENB", "5W59A", "5V5HA",
    "6IF3A", "5DMRA", "4BVXB", "5E8CB", "5OENA", "6LKIB", "6MBAA", "6MAVA",
    "6PNPB", "6ZDTB", "7AZBA", "6FFAA", "3WSOB", "6HLVA", "5JQPB", "4ZV0A",
    "6VCLB", "4ONMB", "4RLJA", "4XAXB", "6UX8A", "5K22B", "5E0QB", "5SZHA",
    "4UYQA", "4W6XA", "5DMBA", "5M2OA", "5V03R", "6DEXA", "6G4JB", "4RLJB",
    "6BCBA", "5C2UA", "6ZIEA", "4Y2OB", "5VGBA", "5NRMA", "6LBUB", "5YNMB",
    "5ZWLE", "6KMQB", "4U5YD", "4ZGYB", "5F22B", "3X37B", "5D1MB", "5IWBA",
    "5MAWE", "5XLYB", "5IWBB", "5LDAA", "6MIBB", "4CMMA", "4CMMB", "5IMMA",
    "4RHZB", "5B64B", "4N80B", "5XBFB", "4ZQUA", "3X37A", "4QJFB", "5WP3B",
    "5ZZAP", "6SM5A", "4BMPB", "6H1FB", "6KP3B", "6H2UB", "6KMQA", "5J26A",
    "5L8LB", "6C7YB", "6OD1B", "3WDGB", "4RS1A", "5JNOA", "5KNHI", "6R6MA",
    "6AW3A", "5GPYB", "6QTAB", "7BXFC", "4PW9B", "5FQ2B", "5NCMB", "4C4KT",
    "5V03B", "6R6MB", "5JWOB", "7CD8A", "4Y5OA", "7BU5A", "5E95B", "6IXXI",
    "5M2OB", "6O6EE", "6WH1A", "4ZV0B", "5KEVA", "6FUBB", "6NDUB", "6ACIH",
    "6SWTB", "4MRTC", "4YIIA", "5F22A", "4UAFE", "5SVHA", "5VGBB", "4UDMA",
    "5B5WU", "5CZDB", "6ODDA", "4K12B", "5DMBD", "5XEDA", "4CJ0B", "6EA3A",
    "6FC0B", "6FUBA", "7D2SB", "4XLGB", "4ZRLB", "5IMTD", "7CN7C", "5CHLA",
    "6ODDB", "6TRIB", "4YH8B", "5J4SB", "5ABUB", "5WXKB", "6BN1B", "6GR8B",
    "5LDAB", "5EU0B", "6KBRC", "4K12A", "4X33A", "7BU5B", "5AIEA", "5JNOB",
    "6BN1A", "5H65B", "6FC3B",
]


BENCHMARK_IDS = {
    "dset186": DSET_186_IDS,
    "pdbtest315": PDBTEST_315_IDS,
}


def parse_benchmark_entry(entry):
    """Return (entry_id, pdb_id, target_chain) for PDB or PDB+chain benchmark IDs."""
    entry_id = entry.strip().upper()
    if len(entry_id) > 4:
        return entry_id, entry_id[:4], entry_id[4:]
    return entry_id, entry_id, None


def download_benchmark_structures(benchmark_entries, output_dir, skip_existing=True):
    """下载标准测试集的PDB结构文件"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    downloaded = []
    for i, entry in enumerate(benchmark_entries):
        entry_id, pdb_id, target_chain = parse_benchmark_entry(entry)
        pdb_path = output_dir / f"{pdb_id}.pdb"
        if skip_existing and pdb_path.exists():
            downloaded.append(
                {
                    "entry_id": entry_id,
                    "pdb_id": pdb_id,
                    "target_chain": target_chain,
                    "pdb_path": pdb_path,
                }
            )
            continue
        try:
            download_and_extract_assembly_pdb(pdb_id, pdb_path)
            downloaded.append(
                {
                    "entry_id": entry_id,
                    "pdb_id": pdb_id,
                    "target_chain": target_chain,
                    "pdb_path": pdb_path,
                }
            )
            print(f"[{i+1:03d}/{len(benchmark_entries)}] Downloaded {entry_id} ({pdb_id})")
        except Exception as exc:
            print(f"[skip] {entry_id} ({pdb_id}): {exc}")
    return downloaded


_DOTS = None  # 全局缓存，避免重复加载


def _get_dots():
    global _DOTS
    if _DOTS is None:
        dot_file = EXAMPLE_DATA_DIR / "Dot.txt"
        _DOTS = load_dots(dot_file)
    return _DOTS


def _summarize_chains(atoms):
    chain_atom_counts = Counter()
    chain_residue_keys = defaultdict(set)
    for atom in atoms:
        chain_atom_counts[atom.chain_id] += 1
        chain_residue_keys[atom.chain_id].add(atom.residue_key())
    chain_ids = sorted(chain_atom_counts)
    chain_residue_counts = {c: len(chain_residue_keys[c]) for c in chain_ids}
    return chain_ids, chain_atom_counts, chain_residue_counts


def process_benchmark_complex(
    pdb_path,
    delta_sasa_threshold=2.0,
    solvent_radius=1.4,
    target_chain=None,
):
    """
    对单个PDB文件用ΔSASA管线计算界面标签，返回残基行列表。
    与 batch_generate_interface_labels.py 的逻辑完全一致。
    """
    try:
        all_atoms = parse_pdb(pdb_path)
        chain_ids, _, chain_residue_counts = _summarize_chains(all_atoms)

        if len(chain_ids) < 2:
            return []

        if target_chain:
            if target_chain not in chain_ids:
                print(f"[warn] {pdb_path.name}: target chain {target_chain} not found")
                return []
            partner_candidates = [
                chain_id
                for chain_id in chain_ids
                if chain_id != target_chain and chain_residue_counts[chain_id] >= 30
            ]
            if not partner_candidates:
                return []
            partner_chain = partner_candidates[0]
        else:
            target_chain = chain_ids[0]
            partner_chain = chain_ids[1]

        if chain_residue_counts[target_chain] < 30 or chain_residue_counts[partner_chain] < 30:
            return []

        dots = _get_dots()
        thresholds = [delta_sasa_threshold]

        # apo SASA：只算 target chain
        target_atoms = filter_atoms_by_chain(all_atoms, target_chain)
        apo_sasa = compute_residue_sasa_for_subset(target_atoms, dots, solvent_radius)

        # holo SASA：target chain 在复合体中的暴露面积
        holo_atoms = filter_atoms_by_chain(all_atoms, [target_chain, partner_chain])
        calculate_sasa(holo_atoms, dots, solvent_radius)
        holo_sasa = aggregate_residue_sasa(
            filter_atoms_by_chain(holo_atoms, target_chain)
        )

        rows = build_interface_rows(apo_sasa, holo_sasa, thresholds)

        # 统一 label 列名，并附加复合体标识（供后续合并数据集使用）
        label_col = f"label_gt_{delta_sasa_threshold:g}"
        pdb_id = pdb_path.stem.upper()
        for row in rows:
            row["label"] = row.get(label_col, 0)
            row["pdb_id"] = pdb_id
            row["target_chain"] = target_chain
            row["partner_chain"] = partner_chain

        return rows

    except Exception as exc:
        print(f"[warn] {pdb_path.name}: {exc}")
        return []


def _compute_sasa_baseline(rows):
    """
    SASA 阈值基线：以 sasa_apo 排名前 K% 的残基作为预测正例（K = 实际正例比例）。
    无需任何模型，纯数据驱动基线。
    """
    labels = [int(row["label"]) for row in rows]
    sasa_vals = [float(row["sasa_apo"]) for row in rows]
    total = len(labels)
    pos = sum(labels)
    if total == 0 or pos == 0:
        return {}

    # 以正例比例作为截断点
    threshold_idx = int(total * pos / total)
    sorted_sasa = sorted(enumerate(sasa_vals), key=lambda x: x[1], reverse=True)
    pred_positive = {idx for idx, _ in sorted_sasa[:pos]}

    tp = sum(1 for i, label in enumerate(labels) if label == 1 and i in pred_positive)
    fp = sum(1 for i, label in enumerate(labels) if label == 0 and i in pred_positive)
    fn = pos - tp
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {"precision": precision, "recall": recall, "f1": f1, "n_pos": pos, "total": total}


def evaluate_on_benchmark(rows, benchmark_name, model_checkpoint=None):
    """
    基准评测：
    - 无 checkpoint：用 SASA 阈值基线估算性能
    - 有 checkpoint：运行模型推断（需数据集已含 ESM 嵌入）
    """
    if not rows:
        print("No rows to evaluate.")
        return

    labels = [int(row["label"]) for row in rows]
    total = len(labels)
    pos = sum(labels)
    neg = total - pos
    print(f"\n=== 基准测试集统计 ===")
    print(f"总残基数: {total}")
    print(f"正例（界面残基）: {pos} ({pos/total*100:.1f}%)")
    print(f"负例（非界面残基）: {neg} ({neg/total*100:.1f}%)")

    # SASA 基线（无需模型）
    baseline = _compute_sasa_baseline(rows)
    if baseline:
        print(f"\n=== SASA 阈值基线（无模型）===")
        print(f"  Precision={baseline['precision']:.4f}  Recall={baseline['recall']:.4f}  F1={baseline['f1']:.4f}")

    if model_checkpoint:
        print(f"\n[checkpoint 推断] 需先提取 ESM-2 嵌入并构建多模态数据集，然后运行：")
        print(f"  python -m sasa_project.train_interface_model \\")
        print(f"    --input data/processed/benchmark_{benchmark_name}_multimodal_650m.csv \\")
        print(f"    --model cross_egnn --predict-only \\")
        print(f"    --manifest data/processed/benchmark_{benchmark_name}_manifest.csv \\")
        print(f"    --checkpoint-input {model_checkpoint}")
    else:
        print("\n=== 完整评测流程（需 GPU 跑 ESM-2）===")
        print("  Step 1 — 提取 ESM-2 嵌入（约 10-30 分钟，需 GPU）:")
        print("    python -m sasa_project.extract_esm_embeddings \\")
        print(f"      --manifest data/processed/benchmark_{benchmark_name}_manifest.csv \\")
        print(f"      --model-name facebook/esm2_t33_650M_UR50D \\")
        print(f"      --device cuda \\")
        print(f"      --output data/processed/benchmark_{benchmark_name}_esm_embeddings_650m.csv")
        print("  Step 2 — 构建多模态数据集:")
        print("    python -m sasa_project.build_multimodal_dataset \\")
        print(f"      --labels data/processed/benchmark_{benchmark_name}_labels.csv \\")
        print(f"      --manifest data/processed/benchmark_{benchmark_name}_manifest.csv \\")
        print(f"      --embeddings data/processed/benchmark_{benchmark_name}_esm_embeddings_650m.csv \\")
        print(f"      --output data/processed/benchmark_{benchmark_name}_multimodal_650m.csv")
        print("  Step 3 — 模型推断:")
        print("    python -m sasa_project.train_interface_model \\")
        print(f"      --input data/processed/benchmark_{benchmark_name}_multimodal_650m.csv \\")
        print("      --model cross_egnn --predict-only \\")
        print(f"      --manifest data/processed/benchmark_{benchmark_name}_manifest.csv \\")
        print("      --checkpoint-input data/processed/best_cross_egnn_650m_esm_sasa_struct_d8.pt")


def save_benchmark_labels(rows, output_path):
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"标签文件已保存: {output_path} ({len(rows)} 行)")


def save_benchmark_manifest(manifest_rows, output_path):
    if not manifest_rows:
        return
    fieldnames = ["entry_id", "pdb_id", "target_chain", "partner_chain", "pdb_path"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Manifest 已保存: {output_path} ({len(manifest_rows)} 个复合体)")


def main():
    parser = argparse.ArgumentParser(
        description="在标准测试集（Dset_186 / PDBtest_315）上评测 PPI 界面预测模型。"
    )
    parser.add_argument(
        "--benchmark",
        choices=["dset186", "pdbtest315"],
        default="dset186",
        help="选择标准测试集",
    )
    parser.add_argument(
        "--output-dir",
        default=str(PDB_COMPLEX_DIR / "benchmarks"),
        help="基准结构文件存放目录",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=2.0,
        help="Delta-SASA 界面定义阈值（A^2）",
    )
    parser.add_argument(
        "--model-checkpoint",
        default=None,
        help="已训练模型的 checkpoint 路径（可选）",
    )
    args = parser.parse_args()

    pdb_ids = BENCHMARK_IDS[args.benchmark]
    print(f"基准测试集: {args.benchmark} ({len(pdb_ids)} 个复合体)")

    # Step 1: 下载结构
    benchmark_items = download_benchmark_structures(pdb_ids, args.output_dir)
    print(f"成功获取 {len(benchmark_items)} 个结构文件")

    # Step 2: 计算 ΔSASA 标签
    all_rows = []
    manifest_rows = []
    for item in benchmark_items:
        pdb_path = item["pdb_path"]
        rows = process_benchmark_complex(
            pdb_path,
            delta_sasa_threshold=args.threshold,
            target_chain=item["target_chain"],
        )
        all_rows.extend(rows)
        if rows:
            first_row = rows[0]
            manifest_rows.append(
                {
                    "entry_id": item["entry_id"],
                    "pdb_id": first_row["pdb_id"],
                    "target_chain": first_row["target_chain"],
                    "partner_chain": first_row["partner_chain"],
                    "pdb_path": str(pdb_path),
                }
            )

    # Step 3: 保存标签文件和 manifest，供 ESM 提取 / 多模态合并复用
    label_output = PROCESSED_DATA_DIR / f"benchmark_{args.benchmark}_labels.csv"
    manifest_output = PROCESSED_DATA_DIR / f"benchmark_{args.benchmark}_manifest.csv"
    save_benchmark_labels(all_rows, label_output)
    save_benchmark_manifest(manifest_rows, manifest_output)

    # Step 4: 评测
    evaluate_on_benchmark(
        all_rows,
        benchmark_name=args.benchmark,
        model_checkpoint=args.model_checkpoint,
    )


if __name__ == "__main__":
    main()
