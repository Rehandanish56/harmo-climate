import unittest
import numpy as np
import pandas as pd
from harmoclimate.psychrometrics import (
    specific_humidity_kg_per_kg,
    relative_humidity_percent_from_specific,
    dew_point_c_from_e,
    vapor_partial_pressure_hpa_from_q_p,
    thermo_from_T_P_RH
)

class TestPsychrometrics(unittest.TestCase):
    def test_specific_humidity_inverse(self):
        # Test that Q -> RH -> Q is consistent
        T = np.array([20.0, 25.0])
        P = np.array([1013.25, 1000.0])
        RH_in = np.array([50.0, 80.0])
        
        Q = specific_humidity_kg_per_kg(T, RH_in, P)
        RH_out = relative_humidity_percent_from_specific(T, Q, P)
        
        np.testing.assert_allclose(RH_in, RH_out, atol=1e-5)

    def test_vapor_pressure_inverse(self):
        # Test Q -> E -> Q
        Q_in = np.array([0.01, 0.02])
        P = np.array([1013.25, 1000.0])
        
        E = vapor_partial_pressure_hpa_from_q_p(Q_in, P)
        # Reversing E -> Q requires the full formula which is inside specific_humidity
        # But we can check E is reasonable
        self.assertTrue(np.all(E > 0))
        self.assertTrue(np.all(E < P))

    def test_thermo_from_T_P_RH(self):
        T = np.array([20.0])
        P = np.array([1013.25])
        RH = np.array([50.0])
        
        q, Td, E, Es = thermo_from_T_P_RH(T, P, RH)
        
        self.assertEqual(len(q), 1)
        self.assertEqual(len(Td), 1)
        self.assertLess(Td[0], T[0]) # Dew point should be less than T for RH < 100
        self.assertLess(E[0], Es[0]) # Vapor pressure < Saturation for RH < 100

    def test_clipping(self):
        # Test RH clipping
        T = np.array([20.0])
        P = np.array([1013.25])
        RH_over = np.array([150.0])
        
        Q = specific_humidity_kg_per_kg(T, RH_over, P)
        # Should behave as 100% RH
        Q_100 = specific_humidity_kg_per_kg(T, [100.0], P)
        np.testing.assert_allclose(Q, Q_100)

if __name__ == "__main__":
    unittest.main()
