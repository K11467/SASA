import argparse
import csv
from pathlib import Path

from delta_sasa_label import (
    build_interface_rows,
    parse_thresholds,
    summarize_thresholds,
    write_interface_labels,
    write_threshold_summary,
)
from sasa import aggregate_residue_sasa, calculate_sasa, filter_atoms_by_chain, load_dots, parse_pdb


def compute_complex_labels(pdb_path, target_chain, partner_chain, dot_file, solvent_radius, thresholds):
    all_atoms = parse_pdb(pdb_path)
    dots = load_dots(dot_file)

    target_atoms = filter_atoms_by_chain(all_atoms, target_chain)
    holo_atoms = filter_atoms_by_chain(all_atoms, [target_chain, partner_chain])

    calculate_sasa(target_atoms, dots, solvent_radius)
    apo_sasa = aggregate_residue_sasa(target_atoms)

    calculate_sasa(holo_atoms, dots, solvent_radius)
    holo_target_atoms = filter_atoms_by_chain(holo_atoms, target_chain)
    holo_sasa = aggregate_residue_sasa(holo_target_atoms)

    rows = build_interface_rows(apo_sasa, holo_sasa, thresholds)
    threshold_rows = summarize_thresholds(rows, thresholds)
    return rows, threshold_rows


def write_aggregate_labels(rows, thresholds, output_path):
    fieldnames = [
        "pdb_id",
        "target_chain",
        "partner_chain",
        "chain_id",
        "residue_id",
        "insertion_code",
        "residue_name",
        "sasa_apo",
        "sasa_holo",
        "delta_sasa",
    ]
    fieldnames.extend(f"label_gt_{threshold:g}" for threshold in thresholds)

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_aggregate_threshold_stats(rows, output_path):
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "pdb_id",
                "target_chain",
                "partner_chain",
                "threshold",
                "positive_count",
                "negative_count",
                "positive_ratio",
                "negative_ratio",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def write_overall_threshold_stats(rows, output_path):
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "threshold",
                "positive_count",
                "negative_count",
                "positive_ratio",
                "negative_ratio",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Batch-generate delta SASA interface labels for a complex dataset."
    )
    parser.add_argument(
        "--manifest",
        default="dataset/complex_manifest.csv",
        help="Manifest CSV produced by collect_complex_dataset.py",
    )
    parser.add_argument("--dot-file", default="Dot.txt", help="Sphere dot file path.")
    parser.add_argument(
        "--solvent-radius", type=float, default=1.4, help="Solvent probe radius."
    )
    parser.add_argument(
        "--thresholds",
        default="0.5,1.0,2.0,5.0",
        help="Comma-separated delta SASA thresholds.",
    )
    parser.add_argument(
        "--per-complex-dir",
        default="dataset/interface_labels_per_complex",
        help="Directory for per-complex label CSV files.",
    )
    parser.add_argument(
        "--per-complex-threshold-dir",
        default="dataset/threshold_stats_per_complex",
        help="Directory for per-complex threshold summary CSV files.",
    )
    parser.add_argument(
        "--aggregate-labels",
        default="dataset/interface_labels_all.csv",
        help="Output CSV for all residue labels.",
    )
    parser.add_argument(
        "--aggregate-thresholds",
        default="dataset/threshold_statistics_by_complex.csv",
        help="Output CSV for threshold stats grouped by complex.",
    )
    parser.add_argument(
        "--overall-thresholds",
        default="dataset/threshold_statistics_overall.csv",
        help="Output CSV for threshold stats over the whole dataset.",
    )
    args = parser.parse_args()

    thresholds = parse_thresholds(args.thresholds)
    per_complex_dir = Path(args.per_complex_dir)
    per_complex_threshold_dir = Path(args.per_complex_threshold_dir)
    per_complex_dir.mkdir(parents=True, exist_ok=True)
    per_complex_threshold_dir.mkdir(parents=True, exist_ok=True)

    aggregate_rows = []
    aggregate_threshold_rows = []
    overall_positive = {threshold: 0 for threshold in thresholds}
    overall_negative = {threshold: 0 for threshold in thresholds}

    with open(args.manifest) as f:
        manifest_rows = list(csv.DictReader(f))

    for index, manifest_row in enumerate(manifest_rows, start=1):
        pdb_id = manifest_row["pdb_id"]
        target_chain = manifest_row["target_chain"]
        partner_chain = manifest_row["partner_chain"]
        pdb_path = manifest_row["pdb_path"]

        rows, threshold_rows = compute_complex_labels(
            pdb_path=pdb_path,
            target_chain=target_chain,
            partner_chain=partner_chain,
            dot_file=args.dot_file,
            solvent_radius=args.solvent_radius,
            thresholds=thresholds,
        )

        labels_path = per_complex_dir / f"{pdb_id}_{target_chain}_{partner_chain}_labels.csv"
        threshold_path = (
            per_complex_threshold_dir / f"{pdb_id}_{target_chain}_{partner_chain}_thresholds.csv"
        )
        write_interface_labels(rows, thresholds, labels_path)
        write_threshold_summary(threshold_rows, threshold_path)

        for row in rows:
            aggregate_row = {
                "pdb_id": pdb_id,
                "target_chain": target_chain,
                "partner_chain": partner_chain,
                "chain_id": row["chain_id"],
                "residue_id": row["residue_id"],
                "insertion_code": row["insertion_code"],
                "residue_name": row["residue_name"],
                "sasa_apo": f"{row['sasa_apo']:.6f}",
                "sasa_holo": f"{row['sasa_holo']:.6f}",
                "delta_sasa": f"{row['delta_sasa']:.6f}",
            }
            for threshold in thresholds:
                aggregate_row[f"label_gt_{threshold:g}"] = row[f"label_gt_{threshold:g}"]
            aggregate_rows.append(aggregate_row)

        for threshold_row in threshold_rows:
            aggregate_threshold_rows.append(
                {
                    "pdb_id": pdb_id,
                    "target_chain": target_chain,
                    "partner_chain": partner_chain,
                    "threshold": f"{threshold_row['threshold']:.2f}",
                    "positive_count": threshold_row["positive_count"],
                    "negative_count": threshold_row["negative_count"],
                    "positive_ratio": f"{threshold_row['positive_ratio']:.6f}",
                    "negative_ratio": f"{threshold_row['negative_ratio']:.6f}",
                }
            )
            threshold = threshold_row["threshold"]
            overall_positive[threshold] += threshold_row["positive_count"]
            overall_negative[threshold] += threshold_row["negative_count"]

        print(
            f"[{index:03d}/{len(manifest_rows)}] "
            f"{pdb_id} {target_chain}-{partner_chain} residues={len(rows)}"
        )

    overall_threshold_rows = []
    for threshold in thresholds:
        positive_count = overall_positive[threshold]
        negative_count = overall_negative[threshold]
        total_count = positive_count + negative_count
        overall_threshold_rows.append(
            {
                "threshold": f"{threshold:.2f}",
                "positive_count": positive_count,
                "negative_count": negative_count,
                "positive_ratio": f"{(positive_count / total_count) if total_count else 0.0:.6f}",
                "negative_ratio": f"{(negative_count / total_count) if total_count else 0.0:.6f}",
            }
        )

    write_aggregate_labels(aggregate_rows, thresholds, args.aggregate_labels)
    write_aggregate_threshold_stats(aggregate_threshold_rows, args.aggregate_thresholds)
    write_overall_threshold_stats(overall_threshold_rows, args.overall_thresholds)

    print(f"Wrote aggregate labels to {args.aggregate_labels}")
    print(f"Wrote per-complex threshold stats to {args.aggregate_thresholds}")
    print(f"Wrote overall threshold stats to {args.overall_thresholds}")


if __name__ == "__main__":
    main()
