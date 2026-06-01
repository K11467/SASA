import math
from collections import OrderedDict, defaultdict
from pathlib import Path

from .paths import PDB_COMPLEX_DIR

# ── 氨基酸疏水性标度（Kyte-Doolittle），用于物理化学特征 ─────────────────────
AA_HYDROPHOBICITY = {
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5,
    "GLN": -3.5, "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5,
    "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8, "PRO": -1.6,
    "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
}


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
}


def resolve_pdb_path(pdb_path):
    path = Path(pdb_path)
    if path.exists():
        return path

    filename = path.name
    fallback_path = PDB_COMPLEX_DIR / filename
    if fallback_path.exists():
        return fallback_path

    raise FileNotFoundError(f"PDB file not found: {pdb_path}")


def residue_sample_id(pdb_id, target_chain, partner_chain, chain_id, residue_id, insertion_code):
    insertion_code = insertion_code or "-"
    return (
        f"{pdb_id}_{target_chain}_{partner_chain}_"
        f"{chain_id}_{residue_id}_{insertion_code}"
    )


def residue_key(chain_id, residue_id, insertion_code):
    return (chain_id, int(residue_id), insertion_code or "")


def residue_key_from_row(row):
    return residue_key(row["chain_id"], row["residue_id"], row.get("insertion_code", ""))


def extract_chain_residues(atoms, chain_id):
    residues = OrderedDict()

    for atom in atoms:
        if atom.chain_id != chain_id:
            continue

        key = residue_key(atom.chain_id, atom.residue_id, atom.insertion_code)
        if key not in residues:
            residues[key] = {
                "chain_id": atom.chain_id,
                "residue_id": atom.residue_id,
                "insertion_code": atom.insertion_code,
                "residue_name": atom.residue_name,
                "atoms": [],
            }
        residues[key]["atoms"].append(atom)

    return list(residues.values())


def residue_sequence(residues):
    return "".join(AA3_TO_1.get(residue["residue_name"].upper(), "X") for residue in residues)


def residue_coordinate(residue):
    atoms = residue["atoms"]
    ca_atoms = [atom for atom in atoms if atom.atom_name == "CA"]
    selected_atoms = ca_atoms if ca_atoms else atoms

    x = sum(atom.x for atom in selected_atoms) / len(selected_atoms)
    y = sum(atom.y for atom in selected_atoms) / len(selected_atoms)
    z = sum(atom.z for atom in selected_atoms) / len(selected_atoms)
    return x, y, z


def residue_coordinate_map(atoms, chain_id):
    coordinate_map = {}
    for residue in extract_chain_residues(atoms, chain_id):
        key = residue_key(
            residue["chain_id"],
            residue["residue_id"],
            residue["insertion_code"],
        )
        coordinate_map[key] = residue_coordinate(residue)
    return coordinate_map


def euclidean_distance(coord_a, coord_b):
    return math.sqrt(
        (coord_a[0] - coord_b[0]) ** 2
        + (coord_a[1] - coord_b[1]) ** 2
        + (coord_a[2] - coord_b[2]) ** 2
    )


# ── Step 3: 二面角 + 半球暴露度特征 ──────────────────────────────────────────

def _dihedral_angle(p0, p1, p2, p3):
    """计算四点定义的二面角（弧度），使用叉积法。"""
    b0 = (p0[0]-p1[0], p0[1]-p1[1], p0[2]-p1[2])
    b1 = (p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2])
    b2 = (p2[0]-p3[0], p2[1]-p3[1], p2[2]-p3[2])

    def cross(a, b):
        return (
            a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0],
        )

    def dot(a, b):
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

    def norm(a):
        return math.sqrt(a[0]**2 + a[1]**2 + a[2]**2)

    n1 = cross(b0, b1)
    n2 = cross(b1, b2)
    n1_norm = norm(n1)
    n2_norm = norm(n2)
    if n1_norm < 1e-8 or n2_norm < 1e-8:
        return 0.0
    n1 = (n1[0]/n1_norm, n1[1]/n1_norm, n1[2]/n1_norm)
    n2 = (n2[0]/n2_norm, n2[1]/n2_norm, n2[2]/n2_norm)
    b1_unit = norm(b1)
    if b1_unit < 1e-8:
        return 0.0
    b1u = (b1[0]/b1_unit, b1[1]/b1_unit, b1[2]/b1_unit)
    x = dot(n1, n2)
    m1 = cross(n1, b1u)
    y = dot(m1, n2)
    return math.atan2(y, x)


def _get_backbone_atoms(residue):
    """从残基原子列表中提取骨架原子坐标（N, CA, C, O, CB）"""
    backbone = {}
    for atom in residue["atoms"]:
        name = atom.atom_name.strip()
        if name in {"N", "CA", "C", "O", "CB"}:
            backbone[name] = (atom.x, atom.y, atom.z)
    return backbone


def compute_backbone_dihedrals(residues):
    """
    计算每个残基的 phi(φ) 和 psi(ψ) 二面角。
    phi_i  = C(i-1) - N(i) - CA(i) - C(i)
    psi_i  = N(i)   - CA(i) - C(i) - N(i+1)
    末端残基无法计算时设为 0。
    返回: list of (sin_phi, cos_phi, sin_psi, cos_psi)，与残基列表等长。
    """
    bb = [_get_backbone_atoms(r) for r in residues]
    n = len(residues)
    result = []
    for i in range(n):
        sin_phi, cos_phi, sin_psi, cos_psi = 0.0, 1.0, 0.0, 1.0

        if i > 0 and {"C"} <= bb[i-1].keys() and {"N", "CA", "C"} <= bb[i].keys():
            phi = _dihedral_angle(bb[i-1]["C"], bb[i]["N"], bb[i]["CA"], bb[i]["C"])
            sin_phi, cos_phi = math.sin(phi), math.cos(phi)

        if i < n-1 and {"N", "CA", "C"} <= bb[i].keys() and {"N"} <= bb[i+1].keys():
            psi = _dihedral_angle(bb[i]["N"], bb[i]["CA"], bb[i]["C"], bb[i+1]["N"])
            sin_psi, cos_psi = math.sin(psi), math.cos(psi)

        result.append((sin_phi, cos_phi, sin_psi, cos_psi))
    return result


def _sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])


def _dot3(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def _norm3(a):
    return math.sqrt(a[0]**2 + a[1]**2 + a[2]**2)


def compute_hse_alpha(residues, radius=13.0):
    """
    计算半球暴露度（HSE-α）。

    算法（Hamelryck 2005）:
      对每个残基 i：
        1. 以 CA(i) 为球心，radius=13Å 为半径的球内的所有其他 CA(j)
        2. 参考向量 v_ref = CA(i) - CA(i-1)（若无前驱则用 CA(i+1) - CA(i)）
        3. 上半球 (HSEup): dot(CA(j)-CA(i), v_ref) > 0 的残基数
        4. 下半球 (HSEdn): dot(CA(j)-CA(i), v_ref) <= 0 的残基数

    返回: list of (hse_up, hse_dn)，与残基列表等长。
    """
    ca_coords = []
    for r in residues:
        ca_list = [a for a in r["atoms"] if a.atom_name.strip() == "CA"]
        if ca_list:
            ca_coords.append((ca_list[0].x, ca_list[0].y, ca_list[0].z))
        else:
            ca_coords.append(None)

    n = len(residues)
    result = []
    for i in range(n):
        if ca_coords[i] is None:
            result.append((0, 0))
            continue

        # 参考向量
        if i > 0 and ca_coords[i-1] is not None:
            v_ref = _sub(ca_coords[i], ca_coords[i-1])
        elif i < n-1 and ca_coords[i+1] is not None:
            v_ref = _sub(ca_coords[i+1], ca_coords[i])
        else:
            result.append((0, 0))
            continue

        ref_len = _norm3(v_ref)
        if ref_len < 1e-8:
            result.append((0, 0))
            continue
        v_ref = (v_ref[0]/ref_len, v_ref[1]/ref_len, v_ref[2]/ref_len)

        hse_up, hse_dn = 0, 0
        for j in range(n):
            if j == i or ca_coords[j] is None:
                continue
            diff = _sub(ca_coords[j], ca_coords[i])
            dist = _norm3(diff)
            if dist > radius:
                continue
            projection = _dot3(diff, v_ref)
            if projection > 0:
                hse_up += 1
            else:
                hse_dn += 1

        result.append((hse_up, hse_dn))
    return result


def extract_structural_features(atoms, chain_id):
    """
    提取指定链的结构特征：二面角（sin/cos）+ 半球暴露度。
    返回: dict，key 为残基 key，value 为特征字典。
    """
    residues = extract_chain_residues(atoms, chain_id)
    dihedrals = compute_backbone_dihedrals(residues)
    hse = compute_hse_alpha(residues)

    features = {}
    for residue, (sin_phi, cos_phi, sin_psi, cos_psi), (hse_up, hse_dn) in zip(
        residues, dihedrals, hse
    ):
        key = residue_key(residue["chain_id"], residue["residue_id"], residue["insertion_code"])
        hydrophobicity = AA_HYDROPHOBICITY.get(residue["residue_name"].upper(), 0.0)
        features[key] = {
            "sin_phi": sin_phi,
            "cos_phi": cos_phi,
            "sin_psi": sin_psi,
            "cos_psi": cos_psi,
            "hse_up": float(hse_up),
            "hse_dn": float(hse_dn),
            "hydrophobicity": hydrophobicity,
        }
    return features


# ─────────────────────────────────────────────────────────────────────────────

def build_edge_pairs(rows, distance_cutoff):
    by_complex = defaultdict(list)
    for index, row in enumerate(rows):
        graph_key = (row["pdb_id"], row["target_chain"], row["partner_chain"])
        by_complex[graph_key].append((index, row))

    edge_pairs = []
    for graph_rows in by_complex.values():
        for left_pos in range(len(graph_rows)):
            left_index, left_row = graph_rows[left_pos]
            left_coord = (
                float(left_row["coord_x"]),
                float(left_row["coord_y"]),
                float(left_row["coord_z"]),
            )
            for right_pos in range(left_pos + 1, len(graph_rows)):
                right_index, right_row = graph_rows[right_pos]
                right_coord = (
                    float(right_row["coord_x"]),
                    float(right_row["coord_y"]),
                    float(right_row["coord_z"]),
                )
                if euclidean_distance(left_coord, right_coord) <= distance_cutoff:
                    edge_pairs.append((left_index, right_index))
                    edge_pairs.append((right_index, left_index))

    return edge_pairs
