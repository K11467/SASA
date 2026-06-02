import argparse
import csv
from pathlib import Path

from .io_utils import read_csv_dicts
from .paths import PROCESSED_DATA_DIR
from .residue_features import (
    extract_structural_features,
    residue_coordinate_map,
    residue_key_from_row,
    residue_sample_id,
    resolve_pdb_path,
)
from .sasa import parse_pdb


def load_embeddings(path):
    if not Path(path).exists():
        return {}

    rows = read_csv_dicts(path)
    return {row["sample_id"]: row for row in rows}


def load_manifest(path):
    rows = read_csv_dicts(path)
    return {
        (row["pdb_id"], row["target_chain"], row["partner_chain"]): row
        for row in rows
    }


def main():
    parser = argparse.ArgumentParser(
        description="Merge delta-SASA weak labels, SASA features, coordinates, and ESM embeddings."
    )
    parser.add_argument(
        "--labels",
        default=str(PROCESSED_DATA_DIR / "interface_labels_all.csv"),
        help="Aggregate delta-SASA label table.",
    )
    parser.add_argument(
        "--embeddings",
        default=str(PROCESSED_DATA_DIR / "esm_residue_embeddings.csv"),
        help="Per-residue ESM embedding CSV.",
    )
    parser.add_argument(
        "--manifest",
        default=str(PROCESSED_DATA_DIR / "complex_manifest.csv"),
        help="Complex manifest CSV.",
    )
    parser.add_argument(
        "--output",
        default=str(PROCESSED_DATA_DIR / "multimodal_residue_dataset.csv"),
        help="Output multimodal residue dataset CSV.",
    )
    parser.add_argument(
        "--label-threshold",
        type=float,
        default=2.0,
        help="Delta-SASA threshold column used as the default label.",
    )
    parser.add_argument(
        "--allow-missing-embeddings",
        action="store_true",
        help="Keep rows without ESM embeddings by filling ESM features with zeros.",
    )
    args = parser.parse_args()

    label_column = f"label_gt_{args.label_threshold:g}"

    label_rows = read_csv_dicts(args.labels)
    if not label_rows:
        raise ValueError(f"Label table is empty: {args.labels}")
    if label_column not in label_rows[0]:
        raise ValueError(f"Label column not found: {label_column}")

    if not Path(args.embeddings).exists() and not args.allow_missing_embeddings:
        raise FileNotFoundError(
            f"Embedding file not found: {args.embeddings}. "
            "Run ESM extraction first or pass --allow-missing-embeddings for SASA-only baselines."
        )

    embeddings_by_id = load_embeddings(args.embeddings)
    embedding_columns = [
        column
        for column in next(iter(embeddings_by_id.values()), {}).keys()
        if column.startswith("esm_")
    ]
    manifest_by_key = load_manifest(args.manifest)
    coordinate_cache = {}
    structural_feature_cache = {}

    output_rows = []
    skipped_missing_embeddings = 0

    for row in label_rows:
        sample_id = residue_sample_id(
            pdb_id=row["pdb_id"],
            target_chain=row["target_chain"],
            partner_chain=row["partner_chain"],
            chain_id=row["chain_id"],
            residue_id=row["residue_id"],
            insertion_code=row.get("insertion_code", ""),
        )
        embedding_row = embeddings_by_id.get(sample_id)
        if embedding_row is None and not args.allow_missing_embeddings:
            skipped_missing_embeddings += 1
            continue

        complex_key = (row["pdb_id"], row["target_chain"], row["partner_chain"])
        if complex_key not in coordinate_cache:
            manifest_row = manifest_by_key[complex_key]
            pdb_path = resolve_pdb_path(manifest_row["pdb_path"])
            atoms = parse_pdb(pdb_path)
            coordinate_cache[complex_key] = residue_coordinate_map(atoms, row["target_chain"])
            structural_feature_cache[complex_key] = extract_structural_features(
                atoms,
                row["target_chain"],
            )

        residue_key = residue_key_from_row(row)
        coord = coordinate_cache[complex_key].get(residue_key)
        if coord is None:
            raise ValueError(f"Missing residue coordinate for sample {sample_id}")
        structural_features = structural_feature_cache[complex_key].get(residue_key, {})

        output_row = {
            "sample_id": sample_id,
            "pdb_id": row["pdb_id"],
            "target_chain": row["target_chain"],
            "partner_chain": row["partner_chain"],
            "chain_id": row["chain_id"],
            "residue_id": row["residue_id"],
            "insertion_code": row["insertion_code"],
            "residue_name": row["residue_name"],
            "coord_x": f"{coord[0]:.6f}",
            "coord_y": f"{coord[1]:.6f}",
            "coord_z": f"{coord[2]:.6f}",
            "sasa_apo": row["sasa_apo"],
            "sasa_holo": row["sasa_holo"],
            "delta_sasa": row["delta_sasa"],
            "label": row[label_column],
            "label_threshold": f"{args.label_threshold:g}",
            "sin_phi": f"{structural_features.get('sin_phi', 0.0):.8f}",
            "cos_phi": f"{structural_features.get('cos_phi', 1.0):.8f}",
            "sin_psi": f"{structural_features.get('sin_psi', 0.0):.8f}",
            "cos_psi": f"{structural_features.get('cos_psi', 1.0):.8f}",
            "hse_up": f"{structural_features.get('hse_up', 0.0):.1f}",
            "hse_dn": f"{structural_features.get('hse_dn', 0.0):.1f}",
            "hydrophobicity": f"{structural_features.get('hydrophobicity', 0.0):.4f}",
        }

        for column in embedding_columns:
            output_row[column] = embedding_row[column] if embedding_row else "0.0"

        output_rows.append(output_row)

    fieldnames = [
        "sample_id",
        "pdb_id",
        "target_chain",
        "partner_chain",
        "chain_id",
        "residue_id",
        "insertion_code",
        "residue_name",
        "coord_x",
        "coord_y",
        "coord_z",
        "sasa_apo",
        "sasa_holo",
        "delta_sasa",
        "label",
        "label_threshold",
        "sin_phi",
        "cos_phi",
        "sin_psi",
        "cos_psi",
        "hse_up",
        "hse_dn",
        "hydrophobicity",
    ]
    fieldnames.extend(embedding_columns)

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    if skipped_missing_embeddings:
        print(f"Skipped {skipped_missing_embeddings} rows without ESM embeddings.")
    print(f"Wrote {len(output_rows)} multimodal residue rows to {args.output}")


if __name__ == "__main__":
    main()
