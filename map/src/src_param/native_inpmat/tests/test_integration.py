from __future__ import annotations

import unittest

import numpy as np

from cama_native_inpmat.integration import _compare_outputs


class IntegrationMetricTests(unittest.TestCase):
    def test_output_comparison_reports_linf_l1_and_missing_mask(self) -> None:
        missing = np.float32(1.0e20)
        explicit = np.array([1.0, 2.0, 3.0, missing], dtype=np.float32)
        composed = np.array([1.0, 2.0, 3.3, missing], dtype=np.float32)
        result = _compare_outputs(explicit, composed)
        self.assertEqual(result["valid_values"], 3)
        self.assertAlmostEqual(result["relative_linf"], 0.1, places=6)
        self.assertAlmostEqual(result["relative_l1"], 0.05, places=6)
        self.assertAlmostEqual(result["max_abs"], 0.3, places=6)


if __name__ == "__main__":
    unittest.main()
