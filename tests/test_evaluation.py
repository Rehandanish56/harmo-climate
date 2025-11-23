
import unittest
import numpy as np
from harmoclimate.training import YearlyDesignStats
from harmoclimate.evaluation import accumulate_climatology_maps, evaluate_loyo

class TestEvaluation(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        self.years = [2000, 2001]
        self.n_samples = 100
        self.feature_dim = 3
        self.stats = []
        
        for year in self.years:
            X = np.random.randn(self.n_samples, self.feature_dim)
            y = np.random.randn(self.n_samples)
            S = X.T @ X
            b = X.T @ y
            
            # Generate valid utc_day (1-365) and utc_hour (0-23)
            # Force overlap by using a small range of days/hours
            utc_day_index = np.random.randint(1, 5, size=self.n_samples)
            utc_hour = np.random.randint(0, 5, size=self.n_samples)
            
            self.stats.append(YearlyDesignStats(
                year=year,
                X=X,
                y=y,
                S=S,
                b=b,
                n=self.n_samples,
                utc_day_index=utc_day_index,
                utc_hour=utc_hour,
                params_meta=[]
            ))

    def test_accumulate_climatology_maps_structure(self):
        total_sum, total_count, yearly_sum, yearly_count = accumulate_climatology_maps(self.stats)
        
        self.assertIsInstance(total_sum, dict)
        self.assertIsInstance(total_count, dict)
        self.assertIsInstance(yearly_sum, dict)
        self.assertIsInstance(yearly_count, dict)
        
        # Check keys are (day, hour) tuples
        if total_sum:
            key = next(iter(total_sum))
            self.assertIsInstance(key, tuple)
            self.assertEqual(len(key), 2)
            self.assertIsInstance(key[0], int)
            self.assertIsInstance(key[1], int)

    def test_evaluate_loyo_runs(self):
        report = evaluate_loyo(self.stats, ridge_lambda=0.0, reference_spec={})
        self.assertIsNotNone(report)
        self.assertEqual(len(report.years), len(self.years))
        self.assertFalse(np.isnan(report.global_rmse))

if __name__ == "__main__":
    unittest.main()
