import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import unittest
import pandas as pd
import numpy as np
from harmoclimate.data_ingest import parse_dt_aaaammjjhh, _process_chunk, StationRecord

class TestDataIngest(unittest.TestCase):
    def test_parse_dt_aaaammjjhh(self):
        # Test clean strings
        s = pd.Series(["2021010100", "2021010101"])
        res = parse_dt_aaaammjjhh(s)
        self.assertEqual(res[0], pd.Timestamp("2021-01-01 00:00:00"))
        self.assertEqual(res[1], pd.Timestamp("2021-01-01 01:00:00"))

        # Test float-like strings
        s = pd.Series(["2021010100.0", "2021010101.0"])
        res = parse_dt_aaaammjjhh(s)
        self.assertEqual(res[0], pd.Timestamp("2021-01-01 00:00:00"))
        
        # Test dirty strings
        s = pd.Series([" 2021010100 ", "invalid"])
        res = parse_dt_aaaammjjhh(s)
        self.assertEqual(res[0], pd.Timestamp("2021-01-01 00:00:00"))
        self.assertTrue(pd.isna(res[1]))

    def test_process_chunk(self):
        # Create a dummy dataframe simulating Meteo-France data
        data = {
            "NUM_POSTE": ["12345678", "12345678.0", "87654321"],
            "NOM_USUEL": ["STATION A", "STATION A", "STATION B"],
            "AAAAMMJJHH": ["2021010112", "2021010113", "2021010112"],
            "T": ["10.5", "11.0", "12.0"],
            "U": ["80", "85", "90"],
            "LON": ["2.35", "2.35", "5.0"],
            "LAT": ["48.85", "48.85", "45.0"],
            "ALTI": ["100", "100", "200"],
            "PSTAT": ["1013.25", "1012.0", "1010.0"]
        }
        df = pd.DataFrame(data)
        
        # Test filtering for station "12345678"
        processed_df, records = _process_chunk(df, "12345678")
        
        self.assertIsNotNone(processed_df)
        self.assertEqual(len(processed_df), 2)
        self.assertEqual(len(records), 2)
        
        # Check processed dataframe columns and types
        self.assertIn("DT_UTC", processed_df.columns)
        self.assertTrue(pd.api.types.is_float_dtype(processed_df["T"]))
        
        # Check records
        self.assertIsInstance(records[0], StationRecord)
        self.assertEqual(records[0].station_code, "12345678")
        self.assertEqual(records[0].station_name, "STATION A")
        self.assertAlmostEqual(records[0].lon, 2.35, places=4)
        
        # Test filtering for non-existent station
        processed_df, records = _process_chunk(df, "99999999")
        self.assertIsNone(processed_df)
        self.assertEqual(len(records), 0)

if __name__ == "__main__":
    unittest.main()
