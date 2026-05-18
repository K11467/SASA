from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw"
PROCESSED_DATA_DIR = PROJECT_ROOT / "data" / "processed"
EXAMPLE_DATA_DIR = RAW_DATA_DIR / "examples"
PDB_COMPLEX_DIR = RAW_DATA_DIR / "pdb_complexes"

