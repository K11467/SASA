import argparse
import csv
import random
from collections import defaultdict

from .residue_features import build_edge_pairs


def require_torch():
    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
    except ImportError as exc:
        raise RuntimeError(
            "Model training requires torch. Install PyTorch first, then rerun this script."
        ) from exc
    return torch, nn, F


def feature_columns(rows, feature_set):
    columns = []
    if feature_set in {"sasa", "esm_sasa"}:
        columns.extend(["sasa_apo", "sasa_holo"])
    if feature_set in {"esm", "esm_sasa"}:
        columns.extend(column for column in rows[0] if column.startswith("esm_"))

    if not columns:
        raise ValueError(f"No feature columns selected for feature set {feature_set!r}.")
    return columns


def split_complexes(rows, seed, train_ratio, val_ratio):
    complex_keys = sorted({(row["pdb_id"], row["target_chain"], row["partner_chain"]) for row in rows})
    random.Random(seed).shuffle(complex_keys)

    train_end = max(1, int(len(complex_keys) * train_ratio))
    val_end = max(train_end + 1, int(len(complex_keys) * (train_ratio + val_ratio)))
    val_end = min(val_end, len(complex_keys))

    train_keys = set(complex_keys[:train_end])
    val_keys = set(complex_keys[train_end:val_end])
    test_keys = set(complex_keys[val_end:])
    if not test_keys:
        test_keys = val_keys
    return train_keys, val_keys, test_keys


def row_complex_key(row):
    return (row["pdb_id"], row["target_chain"], row["partner_chain"])


def select_rows(rows, keys):
    return [row for row in rows if row_complex_key(row) in keys]


def normalize_features(train_x, *other_tensors):
    mean = train_x.mean(dim=0, keepdim=True)
    std = train_x.std(dim=0, keepdim=True)
    std = std.masked_fill(std < 1e-6, 1.0)
    normalized = [(tensor - mean) / std for tensor in (train_x, *other_tensors)]
    return normalized


def tensor_from_rows(torch, rows, columns):
    x = torch.tensor(
        [[float(row[column]) for column in columns] for row in rows],
        dtype=torch.float32,
    )
    y = torch.tensor([float(row["label"]) for row in rows], dtype=torch.float32)
    return x, y


def positive_weight(torch, labels):
    positives = labels.sum()
    negatives = labels.numel() - positives
    if positives <= 0:
        return torch.tensor(1.0)
    return negatives / positives


def binary_metrics(labels, probs):
    pairs = list(zip(labels, probs))
    preds = [1 if prob >= 0.5 else 0 for prob in probs]
    tp = sum(1 for label, pred in zip(labels, preds) if label == 1 and pred == 1)
    tn = sum(1 for label, pred in zip(labels, preds) if label == 0 and pred == 0)
    fp = sum(1 for label, pred in zip(labels, preds) if label == 0 and pred == 1)
    fn = sum(1 for label, pred in zip(labels, preds) if label == 1 and pred == 0)

    accuracy = (tp + tn) / len(labels) if labels else 0.0
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    auroc = compute_auroc(pairs)
    auprc = compute_auprc(pairs)
    return {
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "auroc": auroc,
        "auprc": auprc,
    }


def compute_auroc(pairs):
    positives = [prob for label, prob in pairs if label == 1]
    negatives = [prob for label, prob in pairs if label == 0]
    if not positives or not negatives:
        return 0.0

    wins = 0.0
    total = len(positives) * len(negatives)
    for pos in positives:
        for neg in negatives:
            if pos > neg:
                wins += 1.0
            elif pos == neg:
                wins += 0.5
    return wins / total


def compute_auprc(pairs):
    sorted_pairs = sorted(pairs, key=lambda item: item[1], reverse=True)
    total_positive = sum(1 for label, _ in sorted_pairs if label == 1)
    if total_positive == 0:
        return 0.0

    tp = 0
    fp = 0
    last_recall = 0.0
    area = 0.0
    for label, _ in sorted_pairs:
        if label == 1:
            tp += 1
        else:
            fp += 1
        precision = tp / (tp + fp)
        recall = tp / total_positive
        area += precision * (recall - last_recall)
        last_recall = recall
    return area


class MLPModel:
    def __init__(self, torch, nn, input_dim, hidden_dim, dropout):
        class _MLP(nn.Module):
            def __init__(self):
                super().__init__()
                self.net = nn.Sequential(
                    nn.Linear(input_dim, hidden_dim),
                    nn.ReLU(),
                    nn.Dropout(dropout),
                    nn.Linear(hidden_dim, hidden_dim),
                    nn.ReLU(),
                    nn.Dropout(dropout),
                    nn.Linear(hidden_dim, 1),
                )

            def forward(self, x):
                return self.net(x).squeeze(-1)

        self.model = _MLP()


class GCNModel:
    def __init__(self, torch, nn, F, input_dim, hidden_dim, dropout):
        class _GCN(nn.Module):
            def __init__(self):
                super().__init__()
                self.lin1 = nn.Linear(input_dim, hidden_dim)
                self.lin2 = nn.Linear(hidden_dim, hidden_dim)
                self.out = nn.Linear(hidden_dim, 1)
                self.dropout = dropout

            def forward(self, x, adj):
                x = adj @ self.lin1(x)
                x = F.relu(x)
                x = F.dropout(x, p=self.dropout, training=self.training)
                x = adj @ self.lin2(x)
                x = F.relu(x)
                x = F.dropout(x, p=self.dropout, training=self.training)
                return self.out(x).squeeze(-1)

        self.model = _GCN()


def build_graphs(torch, rows, columns, distance_cutoff):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row_complex_key(row)].append(row)

    graphs = []
    for graph_rows in grouped.values():
        x, y = tensor_from_rows(torch, graph_rows, columns)
        edges = build_edge_pairs(graph_rows, distance_cutoff)
        adj = torch.eye(len(graph_rows), dtype=torch.float32)
        for src, dst in edges:
            adj[src, dst] = 1.0
        degree = adj.sum(dim=1)
        degree_inv_sqrt = degree.pow(-0.5)
        degree_inv_sqrt[torch.isinf(degree_inv_sqrt)] = 0.0
        adj = degree_inv_sqrt.unsqueeze(1) * adj * degree_inv_sqrt.unsqueeze(0)
        graphs.append({"x": x, "y": y, "adj": adj})
    return graphs


def normalize_graphs(torch, train_graphs, *graph_sets):
    train_x = torch.cat([graph["x"] for graph in train_graphs], dim=0)
    mean = train_x.mean(dim=0, keepdim=True)
    std = train_x.std(dim=0, keepdim=True)
    std = std.masked_fill(std < 1e-6, 1.0)
    for graphs in (train_graphs, *graph_sets):
        for graph in graphs:
            graph["x"] = (graph["x"] - mean) / std


def evaluate_mlp(torch, model, x, y):
    model.eval()
    with torch.no_grad():
        probs = torch.sigmoid(model(x)).cpu().tolist()
    labels = [int(value) for value in y.cpu().tolist()]
    return binary_metrics(labels, probs)


def evaluate_gcn(torch, model, graphs):
    model.eval()
    labels = []
    probs = []
    with torch.no_grad():
        for graph in graphs:
            graph_probs = torch.sigmoid(model(graph["x"], graph["adj"]))
            probs.extend(graph_probs.cpu().tolist())
            labels.extend(int(value) for value in graph["y"].cpu().tolist())
    return binary_metrics(labels, probs)


def print_metrics(name, metrics):
    metric_text = " ".join(f"{key}={value:.4f}" for key, value in metrics.items())
    print(f"{name}: {metric_text}")


def train_mlp(torch, nn, rows, columns, split_keys, args):
    train_rows = select_rows(rows, split_keys[0])
    val_rows = select_rows(rows, split_keys[1])
    test_rows = select_rows(rows, split_keys[2])

    train_x, train_y = tensor_from_rows(torch, train_rows, columns)
    val_x, val_y = tensor_from_rows(torch, val_rows, columns)
    test_x, test_y = tensor_from_rows(torch, test_rows, columns)
    train_x, val_x, test_x = normalize_features(train_x, val_x, test_x)

    wrapped = MLPModel(torch, nn, train_x.shape[1], args.hidden_dim, args.dropout)
    model = wrapped.model
    criterion = nn.BCEWithLogitsLoss(pos_weight=positive_weight(torch, train_y))
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    for epoch in range(1, args.epochs + 1):
        model.train()
        optimizer.zero_grad()
        loss = criterion(model(train_x), train_y)
        loss.backward()
        optimizer.step()
        if epoch == 1 or epoch % args.log_every == 0:
            print(f"epoch={epoch:03d} train_loss={loss.item():.6f}")

    print_metrics("train", evaluate_mlp(torch, model, train_x, train_y))
    print_metrics("val", evaluate_mlp(torch, model, val_x, val_y))
    print_metrics("test", evaluate_mlp(torch, model, test_x, test_y))


def train_gcn(torch, nn, F, rows, columns, split_keys, args):
    train_graphs = build_graphs(torch, select_rows(rows, split_keys[0]), columns, args.distance_cutoff)
    val_graphs = build_graphs(torch, select_rows(rows, split_keys[1]), columns, args.distance_cutoff)
    test_graphs = build_graphs(torch, select_rows(rows, split_keys[2]), columns, args.distance_cutoff)
    normalize_graphs(torch, train_graphs, val_graphs, test_graphs)

    wrapped = GCNModel(torch, nn, F, len(columns), args.hidden_dim, args.dropout)
    model = wrapped.model
    train_labels = torch.cat([graph["y"] for graph in train_graphs], dim=0)
    criterion = nn.BCEWithLogitsLoss(pos_weight=positive_weight(torch, train_labels))
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        optimizer.zero_grad()
        for graph in train_graphs:
            loss = criterion(model(graph["x"], graph["adj"]), graph["y"])
            loss.backward()
            total_loss += loss.item()
        optimizer.step()
        if epoch == 1 or epoch % args.log_every == 0:
            print(f"epoch={epoch:03d} train_loss={total_loss / len(train_graphs):.6f}")

    print_metrics("train", evaluate_gcn(torch, model, train_graphs))
    print_metrics("val", evaluate_gcn(torch, model, val_graphs))
    print_metrics("test", evaluate_gcn(torch, model, test_graphs))


def main():
    parser = argparse.ArgumentParser(
        description="Train MLP or lightweight GCN baselines for weakly supervised PPI site prediction."
    )
    parser.add_argument(
        "--input",
        default="data/processed/multimodal_residue_dataset.csv",
        help="Multimodal dataset CSV produced by build_multimodal_dataset.py.",
    )
    parser.add_argument("--model", choices=["mlp", "gcn"], default="mlp")
    parser.add_argument("--feature-set", choices=["sasa", "esm", "esm_sasa"], default="esm_sasa")
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--hidden-dim", type=int, default=128)
    parser.add_argument("--dropout", type=float, default=0.2)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--distance-cutoff", type=float, default=8.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--train-ratio", type=float, default=0.7)
    parser.add_argument("--val-ratio", type=float, default=0.15)
    parser.add_argument("--log-every", type=int, default=5)
    args = parser.parse_args()

    torch, nn, F = require_torch()
    random.seed(args.seed)
    torch.manual_seed(args.seed)

    with open(args.input) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"Input dataset is empty: {args.input}")

    columns = feature_columns(rows, args.feature_set)
    split_keys = split_complexes(rows, args.seed, args.train_ratio, args.val_ratio)
    print(
        f"rows={len(rows)} complexes={sum(len(keys) for keys in split_keys)} "
        f"features={len(columns)} model={args.model} feature_set={args.feature_set}"
    )
    print("Note: delta_sasa is used only to create labels and is not included as a feature.")

    if args.model == "mlp":
        train_mlp(torch, nn, rows, columns, split_keys, args)
    else:
        train_gcn(torch, nn, F, rows, columns, split_keys, args)


if __name__ == "__main__":
    main()
