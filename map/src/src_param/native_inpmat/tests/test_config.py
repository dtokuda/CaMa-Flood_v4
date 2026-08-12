from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import yaml

from cama_native_inpmat.cama import read_diminfo
from cama_native_inpmat.config import (
    configured_output_dir,
    configured_validation_reference_dir,
    load_config,
    resolve_input_files,
    configured_validation_forcing,
)


ROOT = Path(__file__).resolve().parents[1]


class ConfigTests(unittest.TestCase):
    def test_configuration_selects_one_generic_input(self) -> None:
        config = load_config(ROOT / "config" / "cmip6_mrro.yaml")
        self.assertNotIn("cmip_root", config)
        self.assertNotIn("models", config)
        self.assertNotIn("runoff_inputs", config)
        self.assertNotIn("runoff_input", config)
        self.assertEqual(config["output"]["dirname"], "MIROC6")
        self.assertEqual(config["version"], 6)
        self.assertEqual(config["grid"]["coverage"], "global")
        self.assertEqual(config["validation_forcing"], {"use_input": True})
        item = config["input"]
        self.assertEqual(item["variable"], "mrro")
        self.assertTrue(Path(item["path"]).is_absolute())
        self.assertNotIn("output_root", config["cama"])
        self.assertEqual(
            configured_output_dir(config),
            Path(config["cama"]["map_dir"]) / "inpmat_comp-remap" / "MIROC6",
        )
        self.assertEqual(
            configured_validation_reference_dir(config),
            configured_output_dir(config) / "validation_reference",
        )

    def test_reference_mapping_is_explicit_and_matches_diminfo(self) -> None:
        config = load_config(ROOT / "config" / "cmip6_mrro.yaml")
        reference = config["reference_mapping"]
        map_dir = Path(config["cama"]["map_dir"])
        self.assertTrue(Path(reference["diminfo_path"]).is_file())
        self.assertTrue(Path(reference["inpmat_path"]).is_file())
        self.assertEqual(
            Path(reference["diminfo_path"]), map_dir / "diminfo_test-1deg.txt"
        )
        self.assertEqual(
            Path(reference["inpmat_path"]), map_dir / "inpmat_test-1deg.bin"
        )
        diminfo = read_diminfo(reference["diminfo_path"])
        self.assertEqual((diminfo.nxin, diminfo.nyin), (360, 180))
        self.assertEqual(
            reference["input_grid"],
            {
                "type": "global_regular_latlon",
                "longitude_order": "west_to_east",
                "latitude_order": "north_to_south",
            },
        )

    def test_local_input_path_is_nonempty(self) -> None:
        config = load_config(ROOT / "config" / "cmip6_mrro.yaml")
        self.assertTrue(resolve_input_files(config["input"]))

    def test_unknown_configuration_key_is_rejected(self) -> None:
        original = yaml.safe_load((ROOT / "config" / "cmip6_mrro.yaml").read_text())
        original["grid"]["unsupported_option"] = True
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.yaml"
            path.write_text(yaml.safe_dump(original))
            with self.assertRaisesRegex(ValueError, "unknown grid keys"):
                load_config(path)

    def test_mapping_generation_config_may_omit_validation_forcing(self) -> None:
        original = yaml.safe_load((ROOT / "config" / "cmip6_mrro.yaml").read_text())
        original.pop("validation_forcing")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.yaml"
            path.write_text(yaml.safe_dump(original))
            config = load_config(path)
        with self.assertRaisesRegex(ValueError, "requires validation_forcing"):
            configured_validation_forcing(config)


if __name__ == "__main__":
    unittest.main()
