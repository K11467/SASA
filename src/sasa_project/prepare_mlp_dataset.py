import argparse
import csv

from .paths import PROCESSED_DATA_DIR


def threshold_to_column_name(threshold):
    return f"label_gt_{threshold:g}"


def build_sample_id(row):
    insertion_code = row["insertion_code"] or "-"
    return (
        f"{row['pdb_id']}_{row['target_chain']}_{row['partner_chain']}_"
        f"{row['chain_id']}_{row['residue_id']}_{insertion_code}"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Prepare an MLP-friendly residue-level dataset from batch delta SASA labels."
    )
    parser.add_argument(
        "--input",
        default=str(PROCESSED_DATA_DIR / "interface_labels_all.csv"),
        help="Aggregate residue-level interface label table.",
    )
    parser.add_argument(
        "--output",
        default=str(PROCESSED_DATA_DIR / "ml_residue_dataset.csv"),
        help="Output CSV for downstream model training.",
    )
    parser.add_argument(
        "--default-threshold",
        type=float,
        default=2.0,
        help="Threshold used for the default binary label column.",
    )
    args = parser.parse_args()

    default_label_column = threshold_to_column_name(args.default_threshold)

    with open(args.input) as f:
        rows = list(csv.DictReader(f))

    if not rows:
        raise ValueError("Input label table is empty.")

    if default_label_column not in rows[0]:
        raise ValueError(
            f"Default threshold column {default_label_column!r} was not found in {args.input}."
        )

    fieldnames = [
        "sample_id",
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
        "label",
        "label_threshold",
        "label_gt_0.5",
        "label_gt_1",
        "label_gt_2",
        "label_gt_5",
    ]

    output_rows = []
    for row in rows:
        output_rows.append(
            {
                "sample_id": build_sample_id(row),
                "pdb_id": row["pdb_id"],
                "target_chain": row["target_chain"],
                "partner_chain": row["partner_chain"],
                "chain_id": row["chain_id"],
                "residue_id": row["residue_id"],
                "insertion_code": row["insertion_code"],
                "residue_name": row["residue_name"],
                "sasa_apo": row["sasa_apo"],
                "sasa_holo": row["sasa_holo"],
                "delta_sasa": row["delta_sasa"],
                "label": row[default_label_column],
                "label_threshold": f"{args.default_threshold:.1f}",
                "label_gt_0.5": row["label_gt_0.5"],
                "label_gt_1": row["label_gt_1"],
                "label_gt_2": row["label_gt_2"],
                "label_gt_5": row["label_gt_5"],
            }
        )

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    positive_count = sum(int(row["label"]) for row in output_rows)
    total_count = len(output_rows)
    print(
        f"Wrote {total_count} residue samples to {args.output}. "
        f"Default threshold={args.default_threshold:.1f}, "
        f"positive={positive_count}, negative={total_count - positive_count}."
    )


if __name__ == "__main__":
    main()
