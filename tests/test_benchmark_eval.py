import contextlib
import importlib.util
import io
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_benchmark_eval.py"


def load_benchmark_module():
    spec = importlib.util.spec_from_file_location("run_benchmark_eval", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class BenchmarkEvalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_benchmark_module()

    def test_official_pdbtest315_list_is_complete(self):
        self.assertEqual(len(self.module.PDBTEST_315_IDS), 315)

    def test_parse_chain_aware_entry(self):
        self.assertEqual(
            self.module.parse_benchmark_entry("4yoca"),
            ("4YOCA", "4YOC", "A"),
        )

    def test_reproduction_commands_use_selected_benchmark(self):
        rows = [{"label": "1", "sasa_apo": "10.0"}]
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            self.module.evaluate_on_benchmark(rows, benchmark_name="pdbtest315")
        output = stream.getvalue()
        self.assertIn("benchmark_pdbtest315_manifest.csv", output)
        self.assertNotIn("benchmark_dset186_manifest.csv", output)


if __name__ == "__main__":
    unittest.main()
