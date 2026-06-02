import csv
import tempfile
import unittest
from pathlib import Path

from sasa_project.io_utils import read_csv_dicts
from sasa_project.sasa import (
    aggregate_chain_sasa,
    aggregate_residue_sasa,
    calculate_sasa,
    load_dots,
    parse_pdb,
)
from sasa_project.train_interface_model import binary_metrics, feature_columns


ROOT = Path(__file__).resolve().parents[1]


class CsvCompatibilityTests(unittest.TestCase):
    def test_read_csv_dicts_accepts_gbk(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.csv"
            path.write_text("name,value\n测试,1\n", encoding="gbk")
            self.assertEqual(read_csv_dicts(path), [{"name": "测试", "value": "1"}])


class SasaSmokeTests(unittest.TestCase):
    def test_example_sasa_aggregations_match(self):
        atoms = parse_pdb(ROOT / "data" / "raw" / "examples" / "2iww_H.pdb")
        dots = load_dots(ROOT / "data" / "raw" / "examples" / "Dot.txt")
        total = calculate_sasa(atoms, dots)
        residue_total = sum(aggregate_residue_sasa(atoms).values())
        chain_total = sum(aggregate_chain_sasa(atoms).values())
        self.assertAlmostEqual(total, residue_total, places=6)
        self.assertAlmostEqual(total, chain_total, places=6)


class TrainingUtilityTests(unittest.TestCase):
    def test_feature_columns_excludes_delta_sasa(self):
        rows = [{
            "sasa_apo": "1",
            "sasa_holo": "0.5",
            "delta_sasa": "0.5",
            "esm_0": "0.1",
            "sin_phi": "0",
            "cos_phi": "1",
            "sin_psi": "0",
            "cos_psi": "1",
            "hse_up": "1",
            "hse_dn": "2",
            "hydrophobicity": "0.2",
        }]
        columns = feature_columns(rows, "esm_sasa_struct")
        self.assertEqual(len(columns), 10)
        self.assertNotIn("delta_sasa", columns)

    def test_binary_metrics(self):
        metrics = binary_metrics([0, 1, 1, 0], [0.1, 0.9, 0.8, 0.2])
        self.assertEqual(metrics["f1"], 1.0)
        self.assertEqual(metrics["auroc"], 1.0)
        self.assertEqual(metrics["auprc"], 1.0)


if __name__ == "__main__":
    unittest.main()
