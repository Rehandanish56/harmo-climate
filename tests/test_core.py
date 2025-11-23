import unittest
import pandas as pd
import numpy as np
from harmoclimate.core import compute_solar_time, prepare_dataset

class TestCore(unittest.TestCase):
    def test_compute_solar_time_basic(self):
        # Test with simple scalar-like inputs (converted to Series internally)
        dt = pd.to_datetime(["2023-01-01 12:00:00"], utc=True)
        lon = [0.0]
        
        df = compute_solar_time(dt, lon)
        
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 1)
        # At 0 lon, 12:00 UTC, solar hour should be 12.0
        self.assertAlmostEqual(df["hour_solar"].iloc[0], 12.0, places=2)
        self.assertAlmostEqual(df["delta_utc_solar_h"].iloc[0], 0.0, places=2)

    def test_compute_solar_time_series(self):
        # Test with Series inputs
        dt = pd.Series(pd.to_datetime(["2023-01-01 12:00:00", "2023-01-01 12:00:00"], utc=True))
        lon = pd.Series([0.0, 15.0]) # 15 deg = +1 hour
        
        df = compute_solar_time(dt, lon)
        
        self.assertEqual(len(df), 2)
        self.assertAlmostEqual(df["hour_solar"].iloc[0], 12.0, places=2)
        self.assertAlmostEqual(df["hour_solar"].iloc[1], 13.0, places=2) # 12 UTC + 1h offset

    def test_prepare_dataset_missing_columns(self):
        df = pd.DataFrame({"A": [1]})
        with self.assertRaises(KeyError):
            prepare_dataset(df)

    def test_prepare_dataset_calculation(self):
        # Create a dummy dataset
        data = {
            "DT_UTC": pd.to_datetime(["2023-01-01 12:00:00"], utc=True),
            "LON": [0.0],
            "T": [20.0],
            "RH": [50.0],
            "P": [1013.25]
        }
        df = pd.DataFrame(data)
        
        # Request default columns (includes derived ones)
        result = prepare_dataset(df)
        
        self.assertIn("Q", result.columns)
        self.assertIn("E", result.columns)
        self.assertIn("Td", result.columns)
        self.assertIn("hour_solar", result.columns)
        
        # Check values are reasonable
        self.assertGreater(result["Q"].iloc[0], 0)
        self.assertLess(result["Q"].iloc[0], 1) # Specific humidity is small

    def test_prepare_dataset_subset(self):
        data = {
            "DT_UTC": pd.to_datetime(["2023-01-01 12:00:00"], utc=True),
            "LON": [0.0],
            "T": [20.0],
            "RH": [50.0],
            "P": [1013.25]
        }
        df = pd.DataFrame(data)
        
        # Request only specific columns
        result = prepare_dataset(df, columns=["T", "Q"])
        self.assertEqual(set(result.columns), {"T", "Q"})

if __name__ == "__main__":
    unittest.main()
