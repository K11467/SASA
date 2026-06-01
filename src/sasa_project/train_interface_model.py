import argparse
import copy
import csv
import random
from collections import defaultdict
from pathlib import Path

from .residue_features import (
    build_edge_pairs,
    extract_chain_residues,
    residue_coordinate,
)
from .sasa import parse_pdb


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


def resolve_device(torch, requested_device):
    if requested_device == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"
    return requested_device


def feature_columns(rows, feature_set):
    columns = []
    if feature_set in {"sasa", "esm_sasa", "esm_sasa_struct"}:
        columns.extend(["sasa_apo", "sasa_holo"])
    if feature_set in {"esm", "esm_sasa", "esm_sasa_struct"}:
        columns.extend(column for column in rows[0] if column.startswith("esm_"))
    # Step 3: 结构特征（二面角 + 半球暴露度 + 疏水性）
    if feature_set == "esm_sasa_struct":
        struct_cols = ["sin_phi", "cos_phi", "sin_psi", "cos_psi", "hse_up", "hse_dn", "hydrophobicity"]
        columns.extend(c for c in struct_cols if c in rows[0])

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


def to_device(graphs, device):
    for graph in graphs:
        graph["x"] = graph["x"].to(device)
        graph["y"] = graph["y"].to(device)
        graph["adj"] = graph["adj"].to(device)


def positive_weight(torch, labels):
    positives = labels.sum()
    negatives = labels.numel() - positives
    if positives <= 0:
        return torch.tensor(1.0, device=labels.device)
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


# ── Step 4: 等变图神经网络（EGNN）────────────────────────────────────────────
# 论文: E(n) Equivariant Graph Neural Networks, Satorras et al., ICML 2021
# 核心性质: 对旋转、反射、平移具有等变性，符合蛋白质结构的物理对称性。
# 纯 PyTorch 实现，无需额外依赖。

class EGNNLayer:
    """单层 EGNN，在 GCNModel 的 __init__ 闭包风格中定义。"""
    pass


class EGNNModel:
    def __init__(self, torch, nn, F, input_dim, hidden_dim, dropout, n_layers=3):
        class _EGNNLayer(nn.Module):
            """
            EGNN 核心层。
            消息函数: m_ij = φ_e(h_i, h_j, d_ij²)
            坐标更新: x_i += Σ_j (x_i - x_j) · φ_x(m_ij)   [等变]
            特征更新: h_i  = φ_h(h_i, Σ_j m_ij)             [不变]
            """
            def __init__(self, node_dim):
                super().__init__()
                # φ_e: 消息网络，输入 2*node_dim + 1（距离的平方）
                self.phi_e = nn.Sequential(
                    nn.Linear(2 * node_dim + 1, hidden_dim),
                    nn.SiLU(),
                    nn.Linear(hidden_dim, hidden_dim),
                    nn.SiLU(),
                )
                # φ_x: 坐标权重网络，输出标量用于加权方向向量
                self.phi_x = nn.Sequential(
                    nn.Linear(hidden_dim, hidden_dim),
                    nn.SiLU(),
                    nn.Linear(hidden_dim, 1),
                )
                # φ_h: 特征更新网络
                self.phi_h = nn.Sequential(
                    nn.Linear(node_dim + hidden_dim, hidden_dim),
                    nn.SiLU(),
                    nn.Dropout(dropout),
                    nn.Linear(hidden_dim, node_dim),
                )

            def forward(self, h, x, edge_index):
                # h: (N, node_dim), x: (N, 3)
                # edge_index: (2, E) — (src, dst)
                src, dst = edge_index

                # Fix 2: 归一化坐标差，用单位方向向量 + 距离（Å）代替原始差值
                # 避免10-300Å的原始坐标导致dist_sq数值爆炸
                diff = x[src] - x[dst]                        # (E, 3)
                dist = (diff ** 2).sum(dim=-1, keepdim=True).clamp(min=1e-8).sqrt()  # (E, 1)
                diff_norm = diff / dist                        # 单位方向向量
                dist_norm = dist / 10.0                        # 归一化到~[0,1]范围（8Å截断→0.8）

                # 消息函数：用归一化距离而非距离平方
                m_input = torch.cat([h[src], h[dst], dist_norm], dim=-1)
                m = self.phi_e(m_input)                        # (E, hidden_dim)

                # Fix 3: 坐标权重用tanh约束，防止坐标无限漂移
                coord_weight = torch.tanh(self.phi_x(m))      # (E, 1)，值域[-1,1]
                coord_agg = torch.zeros_like(x)
                coord_agg.scatter_add_(0, dst.unsqueeze(-1).expand_as(diff_norm),
                                       diff_norm * coord_weight)
                x = x + coord_agg

                # 特征聚合 + Fix 1: 显式残差连接
                m_agg = torch.zeros(h.shape[0], m.shape[-1], device=h.device)
                m_agg.scatter_add_(0, dst.unsqueeze(-1).expand_as(m), m)
                h = h + self.phi_h(torch.cat([h, m_agg], dim=-1))  # 残差

                return h, x

        class _EGNN(nn.Module):
            def __init__(self):
                super().__init__()
                self.input_proj = nn.Linear(input_dim, hidden_dim)
                self.layers = nn.ModuleList([
                    _EGNNLayer(hidden_dim) for _ in range(n_layers)
                ])
                self.classifier = nn.Sequential(
                    nn.Linear(hidden_dim, hidden_dim // 2),
                    nn.ReLU(),
                    nn.Dropout(dropout),
                    nn.Linear(hidden_dim // 2, 1),
                )

            def forward(self, x_feat, pos, edge_index):
                # x_feat: (N, input_dim), pos: (N, 3), edge_index: (2, E)
                h = self.input_proj(x_feat)
                for layer in self.layers:
                    h, pos = layer(h, pos, edge_index)
                return self.classifier(h).squeeze(-1)

        self.model = _EGNN()


# ── 跨链注意力 EGNN（Cross-Chain EGNN）────────────────────────────────────────
# 核心创新：在 EGNN 输出之后，用 partner chain Cα 坐标做距离权重软注意力，
# 让每个 target chain 残基能直接"感知"配体链的空间布局。
# 物理直觉：界面残基在几何上紧邻 partner chain，距离软注意力可以精准捕捉这一特征。

class CrossChainEGNNModel:
    def __init__(self, torch, nn, F, input_dim, hidden_dim, dropout, n_layers=3):
        class _EGNNLayer(nn.Module):
            def __init__(self, node_dim):
                super().__init__()
                self.phi_e = nn.Sequential(
                    nn.Linear(2 * node_dim + 1, hidden_dim),
                    nn.SiLU(),
                    nn.Linear(hidden_dim, hidden_dim),
                    nn.SiLU(),
                )
                self.phi_x = nn.Sequential(
                    nn.Linear(hidden_dim, hidden_dim),
                    nn.SiLU(),
                    nn.Linear(hidden_dim, 1),
                )
                self.phi_h = nn.Sequential(
                    nn.Linear(node_dim + hidden_dim, hidden_dim),
                    nn.SiLU(),
                    nn.Dropout(dropout),
                    nn.Linear(hidden_dim, node_dim),
                )

            def forward(self, h, x, edge_index):
                src, dst = edge_index
                diff = x[src] - x[dst]
                dist = (diff ** 2).sum(dim=-1, keepdim=True).clamp(min=1e-8).sqrt()
                diff_norm = diff / dist
                dist_norm = dist / 10.0
                m_input = torch.cat([h[src], h[dst], dist_norm], dim=-1)
                m = self.phi_e(m_input)
                coord_weight = torch.tanh(self.phi_x(m))
                coord_agg = torch.zeros_like(x)
                coord_agg.scatter_add_(0, dst.unsqueeze(-1).expand_as(diff_norm),
                                       diff_norm * coord_weight)
                x = x + coord_agg
                m_agg = torch.zeros(h.shape[0], m.shape[-1], device=h.device)
                m_agg.scatter_add_(0, dst.unsqueeze(-1).expand_as(m), m)
                h = h + self.phi_h(torch.cat([h, m_agg], dim=-1))
                return h, x

        class _CrossChainAttention(nn.Module):
            """
            距离权重跨链软注意力。
            Query: target chain 残基 EGNN 特征 h_i
            Key/Value: partner chain Cα 位置投影 partner_h_j
            Attention: alpha_ij = softmax(-d_ij / sigma)  [d_ij = ||pos_i - pos_j||]
            Context: c_i = sum_j alpha_ij * partner_h_j
            """
            def __init__(self):
                super().__init__()
                self.partner_proj = nn.Sequential(
                    nn.Linear(3, hidden_dim),
                    nn.SiLU(),
                    nn.Linear(hidden_dim, hidden_dim),
                )
                # 可学习温度参数，控制注意力软硬程度
                self.log_sigma = nn.Parameter(torch.tensor(2.0))
                self.combine = nn.Sequential(
                    nn.Linear(2 * hidden_dim, hidden_dim),
                    nn.SiLU(),
                    nn.Dropout(dropout),
                )

            def forward(self, h, pos, partner_pos):
                # h: (N, H), pos: (N, 3), partner_pos: (P, 3)
                sigma = self.log_sigma.exp().clamp(min=0.5, max=20.0)
                diff = pos.unsqueeze(1) - partner_pos.unsqueeze(0)   # (N, P, 3)
                dist = diff.norm(dim=-1)                              # (N, P)
                attn = F.softmax(-dist / sigma, dim=-1)               # (N, P)
                partner_h = self.partner_proj(partner_pos)            # (P, H)
                context = attn @ partner_h                            # (N, H)
                return self.combine(torch.cat([h, context], dim=-1))  # (N, H)

        class _CrossChainEGNN(nn.Module):
            def __init__(self):
                super().__init__()
                self.input_proj = nn.Linear(input_dim, hidden_dim)
                self.egnn_layers = nn.ModuleList([
                    _EGNNLayer(hidden_dim) for _ in range(n_layers)
                ])
                self.cross_attn = _CrossChainAttention()
                self.classifier = nn.Sequential(
                    nn.Linear(hidden_dim, hidden_dim // 2),
                    nn.ReLU(),
                    nn.Dropout(dropout),
                    nn.Linear(hidden_dim // 2, 1),
                )

            def forward(self, x_feat, pos, edge_index, partner_pos):
                h = self.input_proj(x_feat)
                for layer in self.egnn_layers:
                    h, pos = layer(h, pos, edge_index)
                h = self.cross_attn(h, pos, partner_pos)
                return self.classifier(h).squeeze(-1)

        self.model = _CrossChainEGNN()


def _load_partner_positions(torch, complex_key, manifest_by_key):
    """从 manifest 加载 partner chain Cα 坐标，失败时返回 None。"""
    try:
        manifest_row = manifest_by_key.get(complex_key)
        if manifest_row is None:
            return None
        from .residue_features import resolve_pdb_path
        pdb_path = resolve_pdb_path(manifest_row["pdb_path"])
        atoms = parse_pdb(str(pdb_path))
        partner_chain = complex_key[2]
        partner_residues = extract_chain_residues(atoms, partner_chain)
        if not partner_residues:
            return None
        coords = [residue_coordinate(r) for r in partner_residues]
        return torch.tensor(coords, dtype=torch.float32)
    except Exception:
        return None


def build_cross_chain_graphs(torch, rows, columns, distance_cutoff, manifest_path):
    """
    与 build_egnn_graphs 相同，但额外加载 partner chain Cα 坐标作为跨链上下文。
    manifest_path: complex_manifest.csv 路径，用于查找 PDB 文件。
    """
    manifest_by_key = {}
    if manifest_path and Path(manifest_path).exists():
        with open(manifest_path) as f:
            for row in csv.DictReader(f):
                key = (row["pdb_id"], row["target_chain"], row["partner_chain"])
                manifest_by_key[key] = row

    grouped = defaultdict(list)
    for row in rows:
        grouped[row_complex_key(row)].append(row)

    graphs = []
    missing_partner = 0
    for complex_key, graph_rows in grouped.items():
        x, y = tensor_from_rows(torch, graph_rows, columns)
        pos = torch.tensor(
            [[float(r["coord_x"]), float(r["coord_y"]), float(r["coord_z"])]
             for r in graph_rows],
            dtype=torch.float32,
        )
        edges = build_edge_pairs(graph_rows, distance_cutoff)
        if edges:
            edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
        else:
            edge_index = torch.zeros((2, 0), dtype=torch.long)

        partner_pos = _load_partner_positions(torch, complex_key, manifest_by_key)
        if partner_pos is None:
            missing_partner += 1
            # 退化：用 target chain 中点作为单个 dummy partner 节点
            partner_pos = pos.mean(dim=0, keepdim=True)

        graphs.append({
            "x": x, "y": y, "pos": pos,
            "edge_index": edge_index,
            "partner_pos": partner_pos,
        })

    if missing_partner:
        print(f"[warn] {missing_partner} complexes missing partner chain positions.")
    return graphs


def to_device_cross_chain(graphs, device):
    for graph in graphs:
        graph["x"] = graph["x"].to(device)
        graph["y"] = graph["y"].to(device)
        graph["pos"] = graph["pos"].to(device)
        graph["edge_index"] = graph["edge_index"].to(device)
        graph["partner_pos"] = graph["partner_pos"].to(device)


def normalize_cross_chain_graphs(torch, train_graphs, *graph_sets):
    train_x = torch.cat([graph["x"] for graph in train_graphs], dim=0)
    mean = train_x.mean(dim=0, keepdim=True)
    std = train_x.std(dim=0, keepdim=True)
    std = std.masked_fill(std < 1e-6, 1.0)
    for graphs in (train_graphs, *graph_sets):
        for graph in graphs:
            graph["x"] = (graph["x"] - mean) / std
    return mean, std


def evaluate_cross_chain(torch, model, graphs):
    model.eval()
    labels, probs = [], []
    with torch.no_grad():
        for graph in graphs:
            out = torch.sigmoid(
                model(graph["x"], graph["pos"], graph["edge_index"], graph["partner_pos"])
            )
            probs.extend(out.cpu().tolist())
            labels.extend(int(v) for v in graph["y"].cpu().tolist())
    return binary_metrics(labels, probs)


def train_cross_chain(torch, nn, F, rows, columns, split_keys, args):
    train_graphs = build_cross_chain_graphs(
        torch, select_rows(rows, split_keys[0]), columns, args.distance_cutoff, args.manifest
    )
    val_graphs = build_cross_chain_graphs(
        torch, select_rows(rows, split_keys[1]), columns, args.distance_cutoff, args.manifest
    )
    test_graphs = build_cross_chain_graphs(
        torch, select_rows(rows, split_keys[2]), columns, args.distance_cutoff, args.manifest
    )

    feature_mean, feature_std = normalize_cross_chain_graphs(
        torch, train_graphs, val_graphs, test_graphs
    )
    to_device_cross_chain(train_graphs, args.device)
    to_device_cross_chain(val_graphs, args.device)
    to_device_cross_chain(test_graphs, args.device)

    wrapped = CrossChainEGNNModel(
        torch, nn, F, len(columns), args.hidden_dim, args.dropout,
        n_layers=getattr(args, "egnn_layers", 3)
    )
    model = wrapped.model.to(args.device)

    train_labels = torch.cat([graph["y"] for graph in train_graphs], dim=0)
    criterion = nn.BCEWithLogitsLoss(pos_weight=positive_weight(torch, train_labels))
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    grad_accum_steps = max(1, args.grad_accum_steps)
    best_val_f1, best_epoch, best_state, stale_epochs = -1.0, 0, None, 0

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        optimizer.zero_grad()
        for graph_index, graph in enumerate(train_graphs, start=1):
            out = model(graph["x"], graph["pos"], graph["edge_index"], graph["partner_pos"])
            loss = criterion(out, graph["y"])
            total_loss += loss.item()
            (loss / grad_accum_steps).backward()
            if graph_index % grad_accum_steps == 0 or graph_index == len(train_graphs):
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=args.max_grad_norm)
                optimizer.step()
                optimizer.zero_grad()
        scheduler.step()

        val_metrics = evaluate_cross_chain(torch, model, val_graphs)
        if val_metrics["f1"] > best_val_f1 + args.early_stop_min_delta:
            best_val_f1 = val_metrics["f1"]
            best_epoch = epoch
            best_state = copy.deepcopy(model.state_dict())
            stale_epochs = 0
        else:
            stale_epochs += 1

        if epoch == 1 or epoch % args.log_every == 0:
            print(
                f"epoch={epoch:03d} "
                f"train_loss={total_loss / max(len(train_graphs), 1):.6f} "
                f"val_f1={val_metrics['f1']:.4f}"
            )

        if args.early_stopping and stale_epochs >= args.patience:
            print(
                f"early_stop epoch={epoch:03d} "
                f"best_epoch={best_epoch:03d} best_val_f1={best_val_f1:.4f}"
            )
            break

    if best_state is not None:
        model.load_state_dict(best_state)

    train_metrics = evaluate_cross_chain(torch, model, train_graphs)
    val_metrics = evaluate_cross_chain(torch, model, val_graphs)
    test_metrics = evaluate_cross_chain(torch, model, test_graphs)
    print_metrics("train", train_metrics)
    print_metrics("val", val_metrics)
    print_metrics("test", test_metrics)

    if args.checkpoint_output:
        metrics = {"train": train_metrics, "val": val_metrics, "test": test_metrics}
        save_checkpoint(
            torch, args.checkpoint_output, model, columns, args, metrics,
            feature_mean=feature_mean, feature_std=feature_std,
        )


def build_egnn_graphs(torch, rows, columns, distance_cutoff):
    """
    为 EGNN 构建图，与 build_graphs 类似，但额外存储 3D 坐标和 edge_index 格式。
    edge_index 格式: (2, E) 的张量，而非对称邻接矩阵。
    """
    grouped = defaultdict(list)
    for row in rows:
        grouped[row_complex_key(row)].append(row)

    graphs = []
    for graph_rows in grouped.values():
        x, y = tensor_from_rows(torch, graph_rows, columns)

        # 3D 坐标
        pos = torch.tensor(
            [[float(r["coord_x"]), float(r["coord_y"]), float(r["coord_z"])]
             for r in graph_rows],
            dtype=torch.float32,
        )

        # 构建 edge_index（src, dst 各一行）
        edges = build_edge_pairs(graph_rows, distance_cutoff)
        if edges:
            edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()  # (2, E)
        else:
            edge_index = torch.zeros((2, 0), dtype=torch.long)

        graphs.append({"x": x, "y": y, "pos": pos, "edge_index": edge_index})
    return graphs


def to_device_egnn(graphs, device):
    for graph in graphs:
        graph["x"] = graph["x"].to(device)
        graph["y"] = graph["y"].to(device)
        graph["pos"] = graph["pos"].to(device)
        graph["edge_index"] = graph["edge_index"].to(device)


def normalize_egnn_graphs(torch, train_graphs, *graph_sets):
    """归一化特征（坐标不归一化，保留等变性）"""
    train_x = torch.cat([graph["x"] for graph in train_graphs], dim=0)
    mean = train_x.mean(dim=0, keepdim=True)
    std = train_x.std(dim=0, keepdim=True)
    std = std.masked_fill(std < 1e-6, 1.0)
    for graphs in (train_graphs, *graph_sets):
        for graph in graphs:
            graph["x"] = (graph["x"] - mean) / std
    return mean, std


def apply_egnn_normalization(torch, graphs, mean, std, device):
    mean = torch.tensor(mean, dtype=torch.float32, device=device)
    std = torch.tensor(std, dtype=torch.float32, device=device)
    for graph in graphs:
        graph["x"] = (graph["x"] - mean.cpu()) / std.cpu()


def evaluate_egnn(torch, model, graphs):
    model.eval()
    labels, probs = [], []
    with torch.no_grad():
        for graph in graphs:
            out = torch.sigmoid(model(graph["x"], graph["pos"], graph["edge_index"]))
            probs.extend(out.cpu().tolist())
            labels.extend(int(v) for v in graph["y"].cpu().tolist())
    return binary_metrics(labels, probs)


def train_egnn(torch, nn, F, rows, columns, split_keys, args):
    train_graphs = build_egnn_graphs(torch, select_rows(rows, split_keys[0]), columns, args.distance_cutoff)
    val_graphs   = build_egnn_graphs(torch, select_rows(rows, split_keys[1]), columns, args.distance_cutoff)
    test_graphs  = build_egnn_graphs(torch, select_rows(rows, split_keys[2]), columns, args.distance_cutoff)

    feature_mean, feature_std = normalize_egnn_graphs(torch, train_graphs, val_graphs, test_graphs)
    to_device_egnn(train_graphs, args.device)
    to_device_egnn(val_graphs, args.device)
    to_device_egnn(test_graphs, args.device)

    wrapped = EGNNModel(torch, nn, F, len(columns), args.hidden_dim, args.dropout,
                        n_layers=getattr(args, "egnn_layers", 3))
    model = wrapped.model.to(args.device)

    train_labels = torch.cat([graph["y"] for graph in train_graphs], dim=0)
    criterion = nn.BCEWithLogitsLoss(pos_weight=positive_weight(torch, train_labels))
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    grad_accum_steps = max(1, args.grad_accum_steps)
    best_val_f1 = -1.0
    best_epoch = 0
    best_state = None
    best_metrics = None
    stale_epochs = 0

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        optimizer.zero_grad()
        for graph_index, graph in enumerate(train_graphs, start=1):
            out = model(graph["x"], graph["pos"], graph["edge_index"])
            loss = criterion(out, graph["y"])
            total_loss += loss.item()

            (loss / grad_accum_steps).backward()
            if graph_index % grad_accum_steps == 0 or graph_index == len(train_graphs):
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=args.max_grad_norm)
                optimizer.step()
                optimizer.zero_grad()
        scheduler.step()

        val_metrics = evaluate_egnn(torch, model, val_graphs)
        if val_metrics["f1"] > best_val_f1 + args.early_stop_min_delta:
            best_val_f1 = val_metrics["f1"]
            best_epoch = epoch
            best_state = copy.deepcopy(model.state_dict())
            best_metrics = {"val": val_metrics}
            stale_epochs = 0
        else:
            stale_epochs += 1

        if epoch == 1 or epoch % args.log_every == 0:
            print(
                f"epoch={epoch:03d} "
                f"train_loss={total_loss / max(len(train_graphs),1):.6f} "
                f"val_f1={val_metrics['f1']:.4f}"
            )

        if args.early_stopping and stale_epochs >= args.patience:
            print(
                f"early_stop epoch={epoch:03d} "
                f"best_epoch={best_epoch:03d} best_val_f1={best_val_f1:.4f}"
            )
            break

    if best_state is not None:
        model.load_state_dict(best_state)

    train_metrics = evaluate_egnn(torch, model, train_graphs)
    val_metrics = evaluate_egnn(torch, model, val_graphs)
    test_metrics = evaluate_egnn(torch, model, test_graphs)
    print_metrics("train", train_metrics)
    print_metrics("val", val_metrics)
    print_metrics("test", test_metrics)

    if args.checkpoint_output:
        metrics = best_metrics or {}
        metrics.update({"train": train_metrics, "val": val_metrics, "test": test_metrics})
        save_checkpoint(
            torch,
            args.checkpoint_output,
            model,
            columns,
            args,
            metrics,
            feature_mean=feature_mean,
            feature_std=feature_std,
        )


def predict_egnn(torch, nn, F, rows, args):
    checkpoint = torch.load(args.checkpoint_input, map_location=args.device)
    if checkpoint.get("model") != "egnn":
        raise ValueError("Only EGNN checkpoints are supported by --predict-only currently.")

    columns = checkpoint["feature_columns"]
    missing_columns = [column for column in columns if column not in rows[0]]
    if missing_columns:
        raise ValueError(f"Prediction input is missing feature columns: {missing_columns[:10]}")
    for coord_column in ("coord_x", "coord_y", "coord_z"):
        if coord_column not in rows[0]:
            raise ValueError(f"Prediction input is missing coordinate column: {coord_column}")

    hidden_dim = checkpoint.get("hidden_dim", args.hidden_dim)
    dropout = checkpoint.get("dropout", args.dropout)
    egnn_layers = checkpoint.get("egnn_layers", args.egnn_layers)
    distance_cutoff = checkpoint.get("distance_cutoff", args.distance_cutoff)

    wrapped = EGNNModel(torch, nn, F, len(columns), hidden_dim, dropout, n_layers=egnn_layers)
    model = wrapped.model.to(args.device)
    model.load_state_dict(checkpoint["model_state_dict"])

    graphs = build_egnn_graphs(torch, rows, columns, distance_cutoff)
    if "feature_mean" not in checkpoint or "feature_std" not in checkpoint:
        raise ValueError("Checkpoint does not contain feature normalization statistics.")
    apply_egnn_normalization(
        torch,
        graphs,
        checkpoint["feature_mean"],
        checkpoint["feature_std"],
        args.device,
    )
    to_device_egnn(graphs, args.device)
    print_metrics("predict", evaluate_egnn(torch, model, graphs))


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


def checkpoint_args(args):
    return {
        key: value
        for key, value in vars(args).items()
        if isinstance(value, (str, int, float, bool, type(None)))
    }


def save_checkpoint(torch, path, model, columns, args, metrics, feature_mean=None, feature_std=None):
    payload = {
        "model_state_dict": model.state_dict(),
        "model": args.model,
        "feature_set": args.feature_set,
        "feature_columns": columns,
        "hidden_dim": args.hidden_dim,
        "dropout": args.dropout,
        "egnn_layers": args.egnn_layers,
        "distance_cutoff": args.distance_cutoff,
        "args": checkpoint_args(args),
        "metrics": metrics,
    }
    if feature_mean is not None and feature_std is not None:
        payload["feature_mean"] = feature_mean.detach().cpu().tolist()
        payload["feature_std"] = feature_std.detach().cpu().tolist()
    torch.save(payload, path)
    print(f"Saved checkpoint: {path}")


def train_mlp(torch, nn, rows, columns, split_keys, args):
    train_rows = select_rows(rows, split_keys[0])
    val_rows = select_rows(rows, split_keys[1])
    test_rows = select_rows(rows, split_keys[2])

    train_x, train_y = tensor_from_rows(torch, train_rows, columns)
    val_x, val_y = tensor_from_rows(torch, val_rows, columns)
    test_x, test_y = tensor_from_rows(torch, test_rows, columns)
    train_x, val_x, test_x = normalize_features(train_x, val_x, test_x)
    train_x = train_x.to(args.device)
    train_y = train_y.to(args.device)
    val_x = val_x.to(args.device)
    val_y = val_y.to(args.device)
    test_x = test_x.to(args.device)
    test_y = test_y.to(args.device)

    wrapped = MLPModel(torch, nn, train_x.shape[1], args.hidden_dim, args.dropout)
    model = wrapped.model.to(args.device)
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
    to_device(train_graphs, args.device)
    to_device(val_graphs, args.device)
    to_device(test_graphs, args.device)

    wrapped = GCNModel(torch, nn, F, len(columns), args.hidden_dim, args.dropout)
    model = wrapped.model.to(args.device)
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
    parser.add_argument("--model", choices=["mlp", "gcn", "egnn", "cross_egnn"], default="mlp")
    parser.add_argument(
        "--feature-set",
        choices=["sasa", "esm", "esm_sasa", "esm_sasa_struct"],
        default="esm_sasa",
        help="esm_sasa_struct 需要数据集中包含二面角/HSE特征（Step 3）",
    )
    parser.add_argument(
        "--egnn-layers",
        type=int,
        default=3,
        help="EGNN 卷积层数（仅对 --model egnn 有效）",
    )
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--hidden-dim", type=int, default=128)
    parser.add_argument("--dropout", type=float, default=0.4)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--distance-cutoff", type=float, default=8.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    parser.add_argument("--train-ratio", type=float, default=0.7)
    parser.add_argument("--val-ratio", type=float, default=0.15)
    parser.add_argument("--log-every", type=int, default=5)
    parser.add_argument(
        "--grad-accum-steps",
        type=int,
        default=8,
        help="EGNN only: accumulate gradients over this many complex graphs before optimizer.step().",
    )
    parser.add_argument(
        "--early-stopping",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="EGNN only: monitor validation F1 and restore the best epoch.",
    )
    parser.add_argument(
        "--patience",
        type=int,
        default=10,
        help="EGNN only: stop after this many epochs without validation F1 improvement.",
    )
    parser.add_argument(
        "--early-stop-min-delta",
        type=float,
        default=1e-4,
        help="EGNN only: minimum validation F1 improvement counted by early stopping.",
    )
    parser.add_argument(
        "--max-grad-norm",
        type=float,
        default=1.0,
        help="EGNN only: gradient clipping norm.",
    )
    parser.add_argument(
        "--checkpoint-output",
        default=None,
        help="Optional path for saving the best model checkpoint.",
    )
    parser.add_argument(
        "--checkpoint-input",
        default=None,
        help="Checkpoint path for --predict-only.",
    )
    parser.add_argument(
        "--predict-only",
        action="store_true",
        help="Load --checkpoint-input and evaluate the rows in --input without training.",
    )
    parser.add_argument(
        "--manifest",
        default="data/processed/complex_manifest.csv",
        help="Complex manifest CSV (required for --model cross_egnn to load partner chain PDB coords).",
    )
    args = parser.parse_args()

    torch, nn, F = require_torch()
    args.device = resolve_device(torch, args.device)
    random.seed(args.seed)
    torch.manual_seed(args.seed)
    print(f"Using device: {args.device}")

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

    if args.predict_only:
        if not args.checkpoint_input:
            raise ValueError("--predict-only requires --checkpoint-input.")
        predict_egnn(torch, nn, F, rows, args)
        return

    if args.model == "mlp":
        train_mlp(torch, nn, rows, columns, split_keys, args)
    elif args.model == "gcn":
        train_gcn(torch, nn, F, rows, columns, split_keys, args)
    elif args.model == "egnn":
        if "coord_x" not in rows[0]:
            raise ValueError(
                "EGNN 需要 3D 坐标列（coord_x/y/z）。"
                "请使用 build_multimodal_dataset.py 生成包含坐标的数据集。"
            )
        train_egnn(torch, nn, F, rows, columns, split_keys, args)
    else:  # cross_egnn
        if "coord_x" not in rows[0]:
            raise ValueError(
                "CrossChainEGNN 需要 3D 坐标列（coord_x/y/z）。"
                "请使用 build_multimodal_dataset.py 生成包含坐标的数据集。"
            )
        print(f"Cross-chain manifest: {args.manifest}")
        train_cross_chain(torch, nn, F, rows, columns, split_keys, args)


if __name__ == "__main__":
    main()
