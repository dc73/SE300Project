import tempfile
import unittest
from pathlib import Path

from fea.computation import run_analysis
from fea.io import export_file, load_file, write_conversion_file


class PipelineTest(unittest.TestCase):
    def test_conversion_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "input.txt"
            nodes = [(1, 0.0, 0.0), (2, 1.0, 0.0), (3, 0.0, 1.0)]
            elements = [(1, 1, 2, 3, 3)]
            material = {"E": 68.9e9, "Poisson": 0.33, "Density": 2700.0}

            self.assertTrue(write_conversion_file(
                nodes, elements, [(1, 1, 1)], [(2, 100, 1)], material, 1, input_path))
            coords, parsed_elements, results = run_analysis(input_path)

            self.assertEqual(coords, [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)])
            self.assertEqual(parsed_elements, [[1, 1, 2, 3, 3]])
            self.assertEqual(len(results["Stress X"]), 3)

    def test_json_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "results.json"
            data = {"results": {"Stress X": [0.0, 1.0]}}
            self.assertTrue(export_file(data, path))
            self.assertEqual(load_file(path), data)


if __name__ == "__main__":
    unittest.main()
