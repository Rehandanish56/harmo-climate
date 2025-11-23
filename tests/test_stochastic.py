import unittest
import numpy as np
import pandas as pd
from harmoclimate.core import simulate_stochastic_year

class TestStochastic(unittest.TestCase):
    def setUp(self):
        # Minimal linear model producing constant value 10.0
        self.linear_model = {
            "model": {
                "n_diurnal": 0,
                "coefficients": [10.0],
                "params_layout": [
                    {"name": "c0", "role": "offset", "n_annual": 0, "start": 0, "length": 1}
                ]
            }
        }

    def test_simulate_var1_deterministic(self):
        # VAR(1) with identity matrix and no noise.
        # But initial state is zero. So it stays zero.
        # T_base = 10. Residuals = 0. T = 10.

        stoch_model = {
            "model": {
                "coefficient_matrix": np.eye(3).tolist(),
                "covariance_matrix": np.zeros((3,3)).tolist(),
                "model_order": 1
            }
        }

        df = simulate_stochastic_year(
            self.linear_model,
            self.linear_model,
            self.linear_model,
            stoch_model,
            year=2023,
            samples_per_day=1 # speed up
        )

        self.assertTrue(np.allclose(df["T"], 10.0))
        self.assertTrue(np.allclose(df["Q"], 10.0))
        self.assertTrue(np.allclose(df["P"], 10.0))

    def test_simulate_var2_runs(self):
        # Just check if it runs with VAR(2) shape
        A = np.zeros((6, 3))
        A[0, 0] = 0.5 # Lag 1 effect on P
        A[3, 0] = 0.2 # Lag 2 effect on P

        stoch_model = {
            "model": {
                "coefficient_matrix": A.tolist(),
                "covariance_matrix": np.eye(3).tolist(), # Some noise
                "model_order": 2
            }
        }

        df = simulate_stochastic_year(
            self.linear_model,
            self.linear_model,
            self.linear_model,
            stoch_model,
            year=2023,
            samples_per_day=1
        )

        self.assertFalse(df.empty)
        # Check that residuals are not all zero (due to noise)
        self.assertFalse(np.allclose(df["T"], 10.0))

if __name__ == "__main__":
    unittest.main()
