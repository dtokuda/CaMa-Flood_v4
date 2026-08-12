from __future__ import annotations

import tempfile
import unittest

from cama_native_inpmat.reporting import write_validation_log


class ReportingTests(unittest.TestCase):
    def test_reset_then_append_uses_one_log_file(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            config = {
                "cama": {"map_dir": directory},
                "output": {"dirname": "example"},
            }
            path = write_validation_log(
                config, "generate_native_inpmat.py", {"status": "passed"}, reset=True
            )
            write_validation_log(
                config, "validate_native_inpmat.py", {"metric": 0.0}
            )
            text = path.read_text(encoding="utf-8")
            self.assertEqual(path.name, "validation.log")
            self.assertIn("generate_native_inpmat.py", text)
            self.assertIn("validate_native_inpmat.py", text)
            self.assertEqual(text.count("status"), 1)


if __name__ == "__main__":
    unittest.main()
