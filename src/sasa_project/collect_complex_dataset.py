import argparse
import csv
import gzip
import json
import subprocess
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

from .paths import PDB_COMPLEX_DIR, PROCESSED_DATA_DIR
from .sasa import parse_pdb


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ASSEMBLY_PDB_GZ_URL_TEMPLATE = "https://files.rcsb.org/download/{pdb_id}.pdb1.gz"

GZIP_MAGIC = b"\x1f\x8b"


def run_curl_json(url, method="GET", data=None):
    command = [
        "curl", "-s", "-L",
        "--retry", "3", "--retry-all-errors",
        "--connect-timeout", "15", "--max-time", "60",
        url,
    ]
    if method == "POST":
        command.extend(["-X", "POST", "-H", "Content-Type: application/json"])
        if data is not None:
            command.extend(["--data", json.dumps(data)])

    result = subprocess.run(command, capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


def download_bytes(url, timeout=30, retries=3):
    """Download URL into memory and return raw bytes. Retries on transient errors."""
    last_error = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "SASA-dataset-collector/1.0"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            if not data:
                raise ValueError("empty response body")
            return data
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, ValueError) as exc:
            last_error = exc
            time.sleep(1.5 * attempt)
    raise RuntimeError(f"download failed after {retries} retries: {last_error}")


def download_and_extract_assembly_pdb(pdb_id, output_path):
    """
    Download .pdb1.gz from RCSB in memory, verify gzip magic bytes,
    decompress, and write atomically to output_path.
    Raises RuntimeError (never crashes the caller) on any failure.
    """
    url = ASSEMBLY_PDB_GZ_URL_TEMPLATE.format(pdb_id=pdb_id)

    try:
        raw = download_bytes(url)
    except RuntimeError as exc:
        raise RuntimeError(f"{pdb_id}: download error: {exc}") from exc

    # Verify magic bytes — HTML error pages start with b'<!'
    if not raw.startswith(GZIP_MAGIC):
        preview = raw[:80]
        raise RuntimeError(
            f"{pdb_id}: response is not gzip (first bytes: {preview!r})"
        )

    try:
        content = gzip.decompress(raw)
    except (gzip.BadGzipFile, OSError) as exc:
        raise RuntimeError(f"{pdb_id}: bad gzip data: {exc}") from exc

    # Atomic write: write to .tmp then rename, so no half-written files
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    tmp_path.write_bytes(content)
    tmp_path.replace(output_path)


def _safe_unlink(path):
    """Delete a file, silently ignoring missing-file and permission errors."""
    try:
        path.unlink(missing_ok=True)
    except OSError:
        pass


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
                        "attribute": "rcsb_assembly_info.polymer_entity_instance_count_protein",
                        "operator": "equals",
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
                        "value": 2.5,
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
                "assembly_id",
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
        "--count", type=int, default=1000, help="Number of valid complexes to collect."
    )
    parser.add_argument(
        "--candidate-limit",
        type=int,
        default=8000,
        help="How many RCSB candidates to inspect before filtering locally.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(PDB_COMPLEX_DIR),
        help="Directory for downloaded PDB files.",
    )
    parser.add_argument(
        "--manifest",
        default=str(PROCESSED_DATA_DIR / "complex_manifest.csv"),
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
    failed = 0

    for pdb_id in candidate_ids:
        if len(manifest_rows) >= args.count:
            break

        pdb_path = output_dir / f"{pdb_id}.pdb"

        # --- download (skip if already on disk) ---
        if not pdb_path.exists():
            try:
                download_and_extract_assembly_pdb(pdb_id, pdb_path)
            except RuntimeError as exc:
                failed += 1
                _safe_unlink(pdb_path)
                print(f"[skip] {pdb_id}: {exc}")
                continue

        # --- local quality filters ---
        try:
            atoms = parse_pdb(pdb_path)
        except Exception as exc:
            failed += 1
            _safe_unlink(pdb_path)
            print(f"[skip] {pdb_id}: parse failed: {exc}")
            continue

        chain_ids, chain_atom_counts, chain_residue_counts = summarize_atoms(atoms)

        if len(chain_ids) != 2:
            rejected += 1
            _safe_unlink(pdb_path)
            continue

        target_chain, partner_chain = chain_ids

        if chain_residue_counts[target_chain] < 30 or chain_residue_counts[partner_chain] < 30:
            rejected += 1
            _safe_unlink(pdb_path)
            continue

        manifest_rows.append(
            {
                "pdb_id": pdb_id,
                "assembly_id": "1",
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

    # Always write manifest with whatever was collected
    write_manifest(manifest_rows, manifest_path)
    print(f"Collected {len(manifest_rows)} valid complexes.")
    print(f"Rejected {rejected} candidates during local filtering.")
    print(f"Failed {failed} candidates due to download or parse errors.")
    print(f"Manifest written to {manifest_path}")

    if len(manifest_rows) < args.count:
        print(
            f"[warn] Only collected {len(manifest_rows)}/{args.count} complexes. "
            f"Try increasing --candidate-limit."
        )


if __name__ == "__main__":
    main()
