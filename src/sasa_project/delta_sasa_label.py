import argparse
import csv

from .paths import EXAMPLE_DATA_DIR, PROCESSED_DATA_DIR
from .sasa import (
    aggregate_residue_sasa,
    calculate_sasa,
    filter_atoms_by_chain,
    get_chain_ids,
    load_dots,
    parse_pdb,
)


DEFAULT_THRESHOLDS = [0.5, 1.0, 2.0, 5.0]


def parse_thresholds(raw_value):
    return [float(value.strip()) for value in raw_value.split(",") if value.strip()]


def compute_residue_sasa_for_subset(atoms, dots, solvent_radius):
    calculate_sasa(atoms, dots, solvent_radius)
    return aggregate_residue_sasa(atoms)


def build_interface_rows(apo_sasa, holo_sasa, thresholds):
    rows = []
    residue_keys = sorted(apo_sasa)

    for key in residue_keys:
        chain_id, residue_id, insertion_code, residue_name = key
        sasa_apo = apo_sasa[key]
        sasa_holo = holo_sasa.get(key, 0.0)
        delta_sasa = sasa_apo - sasa_holo

        row = {
            "chain_id": chain_id,
            "residue_id": residue_id,
            "insertion_code": insertion_code,
            "residue_name": residue_name,
            "sasa_apo": sasa_apo,
            "sasa_holo": sasa_holo,
            "delta_sasa": delta_sasa,
        }

        for threshold in thresholds:
            row[f"label_gt_{threshold:g}"] = 1 if delta_sasa > threshold else 0

        rows.append(row)

    return rows


def write_interface_labels(rows, thresholds, output_path):
    fieldnames = [
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
        for row in rows:
            writer.writerow(
                {
                    "chain_id": row["chain_id"],
                    "residue_id": row["residue_id"],
                    "insertion_code": row["insertion_code"],
                    "residue_name": row["residue_name"],
                    "sasa_apo": f"{row['sasa_apo']:.6f}",
                    "sasa_holo": f"{row['sasa_holo']:.6f}",
                    "delta_sasa": f"{row['delta_sasa']:.6f}",
                    **{
                        f"label_gt_{threshold:g}": row[f"label_gt_{threshold:g}"]
                        for threshold in thresholds
                    },
                }
            )


def summarize_thresholds(rows, thresholds):
    summary_rows = []
    total = len(rows)

    for threshold in thresholds:
        positive = sum(row[f"label_gt_{threshold:g}"] for row in rows)
        negative = total - positive
        positive_ratio = positive / total if total else 0.0
        negative_ratio = negative / total if total else 0.0
        summary_rows.append(
            {
                "threshold": threshold,
                "positive_count": positive,
                "negative_count": negative,
                "positive_ratio": positive_ratio,
                "negative_ratio": negative_ratio,
            }
        )

    return summary_rows


def write_threshold_summary(summary_rows, output_path):
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
        for row in summary_rows:
            writer.writerow(
                {
                    "threshold": f"{row['threshold']:.2f}",
                    "positive_count": row["positive_count"],
                    "negative_count": row["negative_count"],
                    "positive_ratio": f"{row['positive_ratio']:.6f}",
                    "negative_ratio": f"{row['negative_ratio']:.6f}",
                }
            )


def main():
    parser = argparse.ArgumentParser(
        description="Generate interface residue labels from apo/holo delta SASA."
    )
    parser.add_argument(
        "--pdb",
        default=str(EXAMPLE_DATA_DIR / "2WWM.pdb"),
        help="Complex PDB file path.",
    )
    parser.add_argument("--target-chain", required=True, help="Target chain, e.g. A.")
    parser.add_argument(
        "--partner-chains",
        default="",
        help="Partner chains separated by commas. Leave empty to use all non-target chains.",
    )
    parser.add_argument(
        "--dot-file",
        default=str(EXAMPLE_DATA_DIR / "Dot.txt"),
        help="Sphere dot file path.",
    )
    parser.add_argument(
        "--solvent-radius", type=float, default=1.4, help="Solvent probe radius."
    )
    parser.add_argument(
        "--thresholds",
        default="0.5,1.0,2.0,5.0",
        help="Comma-separated delta SASA thresholds.",
    )
    parser.add_argument(
        "--labels-output",
        default=str(PROCESSED_DATA_DIR / "examples" / "interface_labels.csv"),
        help="Per-residue delta SASA labels output CSV.",
    )
    parser.add_argument(
        "--summary-output",
        default=str(PROCESSED_DATA_DIR / "examples" / "threshold_statistics.csv"),
        help="Threshold positive/negative summary CSV.",
    )
    args = parser.parse_args()

    thresholds = parse_thresholds(args.thresholds)
    if not thresholds:
        thresholds = DEFAULT_THRESHOLDS

    all_atoms = parse_pdb(args.pdb)
    available_chains = get_chain_ids(all_atoms)
    if args.target_chain not in available_chains:
        raise ValueError(
            f"Target chain {args.target_chain!r} not found. Available chains: {available_chains}"
        )

    if args.partner_chains.strip():
        partner_chains = [value.strip() for value in args.partner_chains.split(",") if value.strip()]
    else:
        partner_chains = [chain_id for chain_id in available_chains if chain_id != args.target_chain]

    if not partner_chains:
        raise ValueError(
            "No partner chain found for holo calculation. "
            "The current PDB only contains the target chain, so delta SASA would be meaningless."
        )

    invalid_partner_chains = [chain_id for chain_id in partner_chains if chain_id not in available_chains]
    if invalid_partner_chains:
        raise ValueError(
            f"Partner chains not found: {invalid_partner_chains}. Available chains: {available_chains}"
        )

    dots = load_dots(args.dot_file)

    target_atoms = filter_atoms_by_chain(all_atoms, args.target_chain)
    holo_atoms = filter_atoms_by_chain(all_atoms, [args.target_chain, *partner_chains])

    apo_sasa = compute_residue_sasa_for_subset(target_atoms, dots, args.solvent_radius)
    calculate_sasa(holo_atoms, dots, args.solvent_radius)
    holo_sasa = aggregate_residue_sasa(
        filter_atoms_by_chain(holo_atoms, args.target_chain)
    )

    rows = build_interface_rows(apo_sasa, holo_sasa, thresholds)
    summary_rows = summarize_thresholds(rows, thresholds)

    write_interface_labels(rows, thresholds, args.labels_output)
    write_threshold_summary(summary_rows, args.summary_output)

    print(f"Available chains: {available_chains}")
    print(f"Target chain: {args.target_chain}")
    print(f"Partner chains: {partner_chains}")
    print(f"Residue count on target chain: {len(rows)}")
    print(f"Wrote labels to {args.labels_output}")
    print(f"Wrote threshold summary to {args.summary_output}")


if __name__ == "__main__":
    main()
