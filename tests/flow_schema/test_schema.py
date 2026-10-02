"""Catch Flow registration regressions without requiring private Flow access."""
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


class FlowSchemaContract(unittest.TestCase):
    def setUp(self):
        self.schema = json.loads((ROOT / "flow/schema/main.json").read_text())
        self.params = {
            key: spec
            for section in self.schema["inputs"]
            for key, spec in section["params"].items()
        }

    def test_current_flow_input_contract(self):
        self.assertIsInstance(self.schema["inputs"], list)
        for section in self.schema["inputs"]:
            self.assertIn("params", section)
            self.assertNotIn("properties", section)
        for key in ("counts", "samplesheet_", "contrast_table", "gene_sets", "gene_universe"):
            self.assertEqual(self.params[key]["type"], "data", key)
        self.assertEqual(self.params["samplesheet_"]["key"], "samplesheet")
        self.assertTrue(self.params["samplesheet"]["allow_custom_columns"])
        self.assertEqual(self.params["analysis_mode"]["default"], "pairwise")
        self.assertEqual(self.params["module_test"]["default"], "none")
        self.assertFalse((ROOT / "flow/schema/advanced-flow-schema.json").exists())

    def test_statistical_outputs_registered(self):
        for process in ("R_DESIGN_DESEQ2", "R_CAMERA_MODULE"):
            types = {x["filetype"] for x in self.schema["outputs"] if x["process"] == process}
            self.assertTrue({"json", "tsv", "rds"}.issubset(types))

    def test_container_default_matches_pipeline(self):
        default = self.params["custom_model_container"]["default"]
        config = (ROOT / "nextflow.config").read_text()
        self.assertRegex(default, r"^ghcr\.io/[a-z0-9._-]+/differential-analysis/advanced-model@sha256:[0-9a-f]{64}$")
        self.assertIn("custom_model_container = '" + default + "'", config)



if __name__ == "__main__":
    unittest.main()
