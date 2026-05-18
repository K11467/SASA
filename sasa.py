import csv
import math
from collections import defaultdict


VDW_RADII = {
    "C": 1.77,
    "N": 1.66,
    "O": 1.50,
    "S": 1.89,
    "H": 1.20,
}


class Atom:
    def __init__(
        self,
        atom_id,
        atom_name,
        residue_name,
        chain_id,
        residue_id,
        insertion_code,
        x,
        y,
        z,
        element,
    ):
        self.atom_id = atom_id
        self.atom_name = atom_name
        self.residue_name = residue_name
        self.chain_id = chain_id
        self.residue_id = residue_id
        self.insertion_code = insertion_code
        self.x = x
        self.y = y
        self.z = z
        self.element = element
        self.r = VDW_RADII.get(element, 1.50)
        self.r_ext = 0.0
        self.sasa = 0.0

    def residue_key(self):
        return (
            self.chain_id,
            self.residue_id,
            self.insertion_code,
            self.residue_name,
        )


def parse_pdb(filepath):
    atoms = []

    with open(filepath, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            try:
                atom_id = int(line[6:11])
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                residue_id = int(line[22:26])
                insertion_code = line[26:27].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            if not atom_name:
                continue

            element = atom_name[0]
            atoms.append(
                Atom(
                    atom_id=atom_id,
                    atom_name=atom_name,
                    residue_name=residue_name,
                    chain_id=chain_id,
                    residue_id=residue_id,
                    insertion_code=insertion_code,
                    x=x,
                    y=y,
                    z=z,
                    element=element,
                )
            )

    return atoms


def load_dots(filepath):
    dots = []

    with open(filepath, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            dots.append((float(parts[0]), float(parts[1]), float(parts[2])))

    return dots


def calculate_sasa(atoms, dots, solvent_radius=1.4):
    if not atoms or not dots:
        return 0.0

    max_r_ext = 0.0
    for atom in atoms:
        atom.r_ext = atom.r + solvent_radius
        atom.sasa = 0.0
        if atom.r_ext > max_r_ext:
            max_r_ext = atom.r_ext

    cell_size = max_r_ext
    grid = defaultdict(list)

    for i, atom in enumerate(atoms):
        cx = int(atom.x // cell_size)
        cy = int(atom.y // cell_size)
        cz = int(atom.z // cell_size)
        grid[(cx, cy, cz)].append((i, atom))

    total_sasa = 0.0
    num_dots = len(dots)
    neighbor_offsets = [
        (dx, dy, dz)
        for dx in (-1, 0, 1)
        for dy in (-1, 0, 1)
        for dz in (-1, 0, 1)
    ]

    for i, atom in enumerate(atoms):
        exposed_count = 0

        for dx, dy, dz in dots:
            px = atom.x + atom.r_ext * dx
            py = atom.y + atom.r_ext * dy
            pz = atom.z + atom.r_ext * dz

            pcx = int(px // cell_size)
            pcy = int(py // cell_size)
            pcz = int(pz // cell_size)

            is_occluded = False

            for ox, oy, oz in neighbor_offsets:
                neighbor_cell = (pcx + ox, pcy + oy, pcz + oz)
                if neighbor_cell not in grid:
                    continue

                for j, neighbor_atom in grid[neighbor_cell]:
                    if i == j:
                        continue

                    dist_sq = (
                        (px - neighbor_atom.x) ** 2
                        + (py - neighbor_atom.y) ** 2
                        + (pz - neighbor_atom.z) ** 2
                    )
                    if dist_sq < neighbor_atom.r_ext ** 2:
                        is_occluded = True
                        break

                if is_occluded:
                    break

            if not is_occluded:
                exposed_count += 1

        area = (exposed_count / num_dots) * 4 * math.pi * (atom.r_ext ** 2)
        atom.sasa = area
        total_sasa += area

    return total_sasa


def aggregate_residue_sasa(atoms):
    residue_sasa = defaultdict(float)

    for atom in atoms:
        residue_sasa[atom.residue_key()] += atom.sasa

    return residue_sasa


def aggregate_chain_sasa(atoms):
    chain_sasa = defaultdict(float)

    for atom in atoms:
        chain_sasa[atom.chain_id] += atom.sasa

    return chain_sasa


def filter_atoms_by_chain(atoms, chain_ids):
    if isinstance(chain_ids, str):
        chain_ids = [chain_ids]

    chain_ids = set(chain_ids)
    return [atom for atom in atoms if atom.chain_id in chain_ids]


def get_chain_ids(atoms):
    return sorted({atom.chain_id for atom in atoms})


def write_atom_sasa_csv(atoms, output_path):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "atom_id",
                "atom_name",
                "residue_name",
                "chain_id",
                "residue_id",
                "insertion_code",
                "element",
                "x",
                "y",
                "z",
                "sasa",
            ]
        )

        for atom in atoms:
            writer.writerow(
                [
                    atom.atom_id,
                    atom.atom_name,
                    atom.residue_name,
                    atom.chain_id,
                    atom.residue_id,
                    atom.insertion_code,
                    atom.element,
                    f"{atom.x:.3f}",
                    f"{atom.y:.3f}",
                    f"{atom.z:.3f}",
                    f"{atom.sasa:.6f}",
                ]
            )


def write_residue_sasa_csv(residue_sasa, output_path):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "chain_id",
                "residue_id",
                "insertion_code",
                "residue_name",
                "sasa",
            ]
        )

        for key in sorted(residue_sasa):
            chain_id, residue_id, insertion_code, residue_name = key
            writer.writerow(
                [
                    chain_id,
                    residue_id,
                    insertion_code,
                    residue_name,
                    f"{residue_sasa[key]:.6f}",
                ]
            )


def write_chain_sasa_csv(chain_sasa, output_path):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["chain_id", "sasa"])

        for chain_id in sorted(chain_sasa):
            writer.writerow([chain_id, f"{chain_sasa[chain_id]:.6f}"])


def main():
    pdb_file = "2iww_H.pdb"
    dot_file = "Dot.txt"
    solvent_radius = 1.4

    print("Parsing input files...")
    atoms = parse_pdb(pdb_file)
    dots = load_dots(dot_file)
    print(f"Loaded {len(atoms)} atoms and {len(dots)} dots.")

    print("Calculating SASA...")
    total_sasa = calculate_sasa(atoms, dots, solvent_radius)
    residue_sasa = aggregate_residue_sasa(atoms)
    chain_sasa = aggregate_chain_sasa(atoms)

    write_atom_sasa_csv(atoms, "atom_sasa.csv")
    write_residue_sasa_csv(residue_sasa, "residue_sasa.csv")
    write_chain_sasa_csv(chain_sasa, "chain_sasa.csv")

    atom_sum = sum(atom.sasa for atom in atoms)
    residue_sum = sum(residue_sasa.values())

    print(f"Total SASA: {total_sasa:.6f}")
    print(f"Atom count: {len(atoms)}")
    print(f"Residue count: {len(residue_sasa)}")
    print(f"Chain count: {len(chain_sasa)}")
    print(f"Check total == sum(atom.sasa): {atom_sum:.6f}")
    print(f"Check total == sum(residue_sasa): {residue_sum:.6f}")
    print("Wrote atom_sasa.csv, residue_sasa.csv, chain_sasa.csv")


if __name__ == "__main__":
    main()
