import argparse
import csv
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from sasa import parse_pdb


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY_URL_TEMPLATE = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
PDB_URL_TEMPLATE = "https://files.rcsb.org/download/{pdb_id}.pdb"


def run_curl_json(url, method="GET", data=None):
    command = ["curl", "-s", "-L", "--retry", "3", "--retry-all-errors", url]
    if method == "POST":
        command.extend(["-X", "POST", "-H", "Content-Type: application/json"])
        if data is not None:
            command.extend(["--data", json.dumps(data)])

    result = subprocess.run(command, capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


def download_file(url, output_path):
    subprocess.run(
        ["curl", "-s", "-L", "--retry", "3", "--retry-all-errors", url, "-o", str(output_path)],
        check=True,
        capture_output=True,
        text=True,
    )


def build_search_query(limit):
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.selected_polymer_entity_types",
                        "operator": "exact_match",
                        "value": "Protein (only)",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        "operator": "greater_or_equal",
                        "value": 2,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": 3.0,
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": limit},
            "results_content_type": ["experimental"],
            "sort": [
                {
                    "sort_by": "rcsb_entry_info.resolution_combined",
                    "direction": "asc",
                }
            ],
        },
    }


def fetch_candidate_ids(limit):
    payload = build_search_query(limit)
    result = run_curl_json(SEARCH_URL, method="POST", data=payload)
    return [item["identifier"] for item in result.get("result_set", [])]


def fetch_entry_metadata(pdb_id):
    metadata = run_curl_json(ENTRY_URL_TEMPLATE.format(pdb_id=pdb_id))
    title = metadata.get("struct", {}).get("title", "")
    resolutions = metadata.get("rcsb_entry_info", {}).get("resolution_combined", [])
    resolution = resolutions[0] if resolutions else ""
    return {
        "title": title,
        "resolution": resolution,
    }


def summarize_atoms(atoms):
    chain_atom_counts = Counter()
    chain_residue_keys = defaultdict(set)

    for atom in atoms:
        chain_atom_counts[atom.chain_id] += 1
        chain_residue_keys[atom.chain_id].add(atom.residue_key())

    chain_ids = sorted(chain_atom_counts)
    chain_residue_counts = {
        chain_id: len(chain_residue_keys[chain_id]) for chain_id in chain_ids
    }
    return chain_ids, chain_atom_counts, chain_residue_counts


def write_manifest(rows, output_path):
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "pdb_id",
                "title",
                "resolution",
                "target_chain",
                "partner_chain",
                "target_residue_count",
                "partner_residue_count",
                "target_atom_count",
                "partner_atom_count",
                "pdb_path",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Collect a clean two-chain protein complex dataset from RCSB PDB."
    )
    parser.add_argument(
        "--count", type=int, default=100, help="Number of valid complexes to collect."
    )
    parser.add_argument(
        "--candidate-limit",
        type=int,
        default=600,
        help="How many RCSB candidates to inspect before filtering locally.",
    )
    parser.add_argument(
        "--output-dir",
        default="dataset/pdb_complexes",
        help="Directory for downloaded PDB files.",
    )
    parser.add_argument(
        "--manifest",
        default="dataset/complex_manifest.csv",
        help="Output CSV for collected complex metadata.",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(args.manifest)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    candidate_ids = fetch_candidate_ids(args.candidate_limit)
    print(f"Fetched {len(candidate_ids)} candidate PDB entries from RCSB.")

    manifest_rows = []
    rejected = 0

    for pdb_id in candidate_ids:
        if len(manifest_rows) >= args.count:
            break

        pdb_path = output_dir / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            download_file(PDB_URL_TEMPLATE.format(pdb_id=pdb_id), pdb_path)

        atoms = parse_pdb(pdb_path)
        chain_ids, chain_atom_counts, chain_residue_counts = summarize_atoms(atoms)

        if len(chain_ids) != 2:
            rejected += 1
            pdb_path.unlink(missing_ok=True)
            continue

        target_chain, partner_chain = chain_ids

        if chain_residue_counts[target_chain] < 20 or chain_residue_counts[partner_chain] < 20:
            rejected += 1
            pdb_path.unlink(missing_ok=True)
            continue

        metadata = fetch_entry_metadata(pdb_id)
        manifest_rows.append(
            {
                "pdb_id": pdb_id,
                "title": metadata["title"],
                "resolution": metadata["resolution"],
                "target_chain": target_chain,
                "partner_chain": partner_chain,
                "target_residue_count": chain_residue_counts[target_chain],
                "partner_residue_count": chain_residue_counts[partner_chain],
                "target_atom_count": chain_atom_counts[target_chain],
                "partner_atom_count": chain_atom_counts[partner_chain],
                "pdb_path": str(pdb_path),
            }
        )
        print(
            f"[{len(manifest_rows):03d}/{args.count}] "
            f"{pdb_id} {target_chain}-{partner_chain} "
            f"residues={chain_residue_counts[target_chain]}/{chain_residue_counts[partner_chain]}"
        )

    write_manifest(manifest_rows, manifest_path)
    print(f"Collected {len(manifest_rows)} valid complexes.")
    print(f"Rejected {rejected} candidates during local filtering.")
    print(f"Manifest written to {manifest_path}")

    if len(manifest_rows) < args.count:
        raise RuntimeError(
            f"Only collected {len(manifest_rows)} complexes, fewer than requested {args.count}."
        )


if __name__ == "__main__":
    main()
