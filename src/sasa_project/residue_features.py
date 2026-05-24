import math
from collections import OrderedDict, defaultdict
from pathlib import Path

from .paths import PDB_COMPLEX_DIR


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
