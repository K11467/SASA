import argparse
import csv

from .io_utils import read_csv_dicts
from .paths import PROCESSED_DATA_DIR
from .residue_features import (
    extract_chain_residues,
    residue_key,
    residue_sample_id,
    residue_sequence,
    resolve_pdb_path,
)
from .sasa import parse_pdb


def resolve_device(torch, requested_device):
    if requested_device != "auto":
        return requested_device
    return "cuda" if torch.cuda.is_available() else "cpu"


def load_model(model_name, device):
    try:
        import torch
        from transformers import AutoModel, AutoTokenizer
    except ImportError as exc:
        raise RuntimeError(
            "ESM embedding extraction requires torch and transformers. "
            "Install them first, then rerun this script."
        ) from exc

    device = resolve_device(torch, device)
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)
    model.eval()
    try:
        model.to(device)
    except Exception as exc:
        if device == "cuda":
            print(f"CUDA device is not usable for this PyTorch build: {exc}")
            print("Falling back to CPU.")
            device = "cpu"
            model.to(device)
        else:
            raise
    return torch, tokenizer, model, device


def extract_embeddings_for_sequence(torch, tokenizer, model, sequence, device):
    inputs = tokenizer(sequence, return_tensors="pt", add_special_tokens=True)
    inputs = {key: value.to(device) for key, value in inputs.items()}

    with torch.no_grad():
        outputs = model(**inputs)

    token_embeddings = outputs.last_hidden_state[0]
    residue_embeddings = token_embeddings[1 : len(sequence) + 1].detach().cpu().tolist()
    if len(residue_embeddings) != len(sequence):
        raise ValueError(
            f"Expected {len(sequence)} residue embeddings, got {len(residue_embeddings)}."
        )
    return residue_embeddings


def main():
    parser = argparse.ArgumentParser(
        description="Extract per-residue ESM-2 embeddings for target chains in the complex manifest."
    )
    parser.add_argument(
        "--manifest",
        default=str(PROCESSED_DATA_DIR / "complex_manifest.csv"),
        help="Complex manifest CSV.",
    )
    parser.add_argument(
        "--output",
        default=str(PROCESSED_DATA_DIR / "esm_residue_embeddings.csv"),
        help="Output CSV containing one ESM embedding vector per residue.",
    )
    parser.add_argument(
        "--model-name",
        default="facebook/esm2_t6_8M_UR50D",
        help="HuggingFace model id or local model path.",
    )
    parser.add_argument(
        "--device",
        default="auto",
        help="Torch device: auto, cpu, or cuda.",
    )
    parser.add_argument(
        "--max-complexes",
        type=int,
        default=0,
        help="Limit processed complexes for quick tests. 0 means all.",
    )
    args = parser.parse_args()

    torch, tokenizer, model, device = load_model(args.model_name, args.device)

    manifest_rows = read_csv_dicts(args.manifest)
    if args.max_complexes:
        manifest_rows = manifest_rows[: args.max_complexes]

    embedding_dim = None
    writer = None
    output_file = open(args.output, "w", newline="")

    try:
        for index, manifest_row in enumerate(manifest_rows, start=1):
            pdb_id = manifest_row["pdb_id"]
            target_chain = manifest_row["target_chain"]
            partner_chain = manifest_row["partner_chain"]
            pdb_path = resolve_pdb_path(manifest_row["pdb_path"])

            atoms = parse_pdb(pdb_path)
            residues = extract_chain_residues(atoms, target_chain)
            sequence = residue_sequence(residues)
            if not sequence:
                raise ValueError(f"No residues found for {pdb_id} chain {target_chain}.")

            embeddings = extract_embeddings_for_sequence(
                torch=torch,
                tokenizer=tokenizer,
                model=model,
                sequence=sequence,
                device=device,
            )
            embedding_dim = len(embeddings[0])

            if writer is None:
                fieldnames = [
                    "sample_id",
                    "pdb_id",
                    "target_chain",
                    "partner_chain",
                    "chain_id",
                    "residue_id",
                    "insertion_code",
                    "residue_name",
                ]
                fieldnames.extend(f"esm_{dim}" for dim in range(embedding_dim))
                writer = csv.DictWriter(output_file, fieldnames=fieldnames)
                writer.writeheader()

            for residue, embedding in zip(residues, embeddings):
                output_row = {
                    "sample_id": residue_sample_id(
                        pdb_id=pdb_id,
                        target_chain=target_chain,
                        partner_chain=partner_chain,
                        chain_id=residue["chain_id"],
                        residue_id=residue["residue_id"],
                        insertion_code=residue["insertion_code"],
                    ),
                    "pdb_id": pdb_id,
                    "target_chain": target_chain,
                    "partner_chain": partner_chain,
                    "chain_id": residue["chain_id"],
                    "residue_id": residue["residue_id"],
                    "insertion_code": residue["insertion_code"],
                    "residue_name": residue["residue_name"],
                }
                for dim, value in enumerate(embedding):
                    output_row[f"esm_{dim}"] = f"{value:.8f}"
                writer.writerow(output_row)

            print(
                f"[{index:03d}/{len(manifest_rows)}] "
                f"{pdb_id} chain={target_chain} residues={len(residues)} "
                f"dim={embedding_dim}"
            )
    finally:
        output_file.close()

    print(f"Wrote residue embeddings to {args.output}")


if __name__ == "__main__":
    main()
