import argparse
import csv
from pathlib import Path

from sasa_project.io_utils import read_csv_dicts
from sasa_project.train_interface_model import binary_metrics, select_rows, split_complexes


def evaluate_apo_rank_baseline(rows):
    labels = [int(row["label"]) for row in rows]
    positive_count = sum(labels)
    ranked = sorted(
        enumerate(rows),
        key=lambda item: float(item[1]["sasa_apo"]),
        reverse=True,
    )
    predicted_positive = {index for index, _ in ranked[:positive_count]}
    probabilities = [1.0 if index in predicted_positive else 0.0 for index in range(len(rows))]
    return binary_metrics(labels, probabilities)


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate an apo-SASA rank baseline on the same complex split as trained models."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--train-ratio", type=float, default=0.7)
    parser.add_argument("--val-ratio", type=float, default=0.15)
    args = parser.parse_args()

    rows = read_csv_dicts(args.input)
    split_keys = split_complexes(rows, args.seed, args.train_ratio, args.val_ratio)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["split", "model", "feature_set", "accuracy", "precision", "recall", "f1", "auroc", "auprc"]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for split, keys in zip(("train", "val", "test"), split_keys):
            metrics = evaluate_apo_rank_baseline(select_rows(rows, keys))
            writer.writerow({
                "split": split,
                "model": "apo_sasa_rank",
                "feature_set": "apo",
                **{key: f"{value:.6f}" for key, value in metrics.items()},
            })
            print(split, " ".join(f"{key}={value:.4f}" for key, value in metrics.items()))


if __name__ == "__main__":
    main()
