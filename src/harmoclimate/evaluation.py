"""Evaluation protocols for HarmoClimate linear models.

Leave-one-year-out (LOYO) validation now evaluates on a no-leap UTC grid while
models remain trained on the solar feature space. The climatology baseline
averages each (utc_day_of_year, utc_hour) bucket excluding the held-out year
and aggregates errors using observation-count weighting.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Dict, List, Tuple

import numpy as np

from .training import (
    FINAL_TRAINING_PERIOD_LABEL,
    LeaveOneYearOutReport,
    YearlyDesignStats,
    YearlyValidationMetrics,
    solve_normal_equations,
)


def _accumulate_climatology_arrays(
    stats: List[YearlyDesignStats],
) -> Tuple[
    np.ndarray,
    np.ndarray,
    Dict[int, Tuple[np.ndarray, np.ndarray]],
]:
    """Accumulate per-grid sums/counts using numpy arrays.
    
    Returns:
        total_sum_arr: (367, 24) float array
        total_count_arr: (367, 24) int array
        yearly_data: Dict[year, (sum_arr, count_arr)]
    """
    # Shape: (367, 24) to cover days 1-366 (index 0 unused) and hours 0-23
    shape = (367, 24)
    total_sum = np.zeros(shape, dtype=float)
    total_count = np.zeros(shape, dtype=int)
    yearly_data = {}

    for stat in stats:
        # Filter valid indices
        valid = (
            (stat.utc_day_index >= 1)
            & (stat.utc_day_index <= 366)
            & (stat.utc_hour >= 0)
            & (stat.utc_hour < 24)
            & np.isfinite(stat.y)
        )
        
        days = stat.utc_day_index[valid]
        hours = stat.utc_hour[valid]
        vals = stat.y[valid]
        
        # Accumulate for this year
        year_sum = np.zeros(shape, dtype=float)
        year_count = np.zeros(shape, dtype=int)
        
        np.add.at(year_sum, (days, hours), vals)
        np.add.at(year_count, (days, hours), 1)
        
        yearly_data[stat.year] = (year_sum, year_count)
        
        # Accumulate to total
        total_sum += year_sum
        total_count += year_count
        
    return total_sum, total_count, yearly_data


def accumulate_climatology_maps(
    stats: List[YearlyDesignStats],
) -> Tuple[
    Dict[Tuple[int, int], float],
    Dict[Tuple[int, int], int],
    Dict[int, Dict[Tuple[int, int], float]],
    Dict[int, Dict[Tuple[int, int], int]],
]:
    """Accumulate per-grid sums/counts for the UTC climatology baseline.

    Returns dictionaries for (utc_day_index, utc_hour) → value, keeping both the
    global aggregates and per-year contributions to support held-out exclusion.
    """
    total_sum_arr, total_count_arr, yearly_data_arr = _accumulate_climatology_arrays(stats)
    
    # Convert arrays to dicts for backward compatibility
    def arr_to_dict(arr, dtype):
        d = {}
        rows, cols = np.nonzero(arr)
        for r, c in zip(rows, cols):
            d[(int(r), int(c))] = dtype(arr[r, c])
        return d

    total_sum = arr_to_dict(total_sum_arr, float)
    total_count = arr_to_dict(total_count_arr, int)
    
    yearly_sum = {}
    yearly_count = {}
    
    for year, (y_sum, y_count) in yearly_data_arr.items():
        yearly_sum[year] = arr_to_dict(y_sum, float)
        yearly_count[year] = arr_to_dict(y_count, int)

    return total_sum, total_count, yearly_sum, yearly_count


def evaluate_loyo(
    stats: List[YearlyDesignStats],
    *,
    ridge_lambda: float,
    reference_spec: Dict[str, object],
) -> LeaveOneYearOutReport:
    """Run leave-one-year-out validation with a UTC climatology baseline.

    The baseline computes a mean per (utc_day_index, utc_hour) cell on a
    no-leap (1–365) grid and excludes the held-out year before each
    comparison. Aggregation of RMSE/skill statistics is weighted by the number
    of valid observations per fold.
    """

    if not stats:
        return LeaveOneYearOutReport(
            years=[],
            global_rmse=math.nan,
            global_skill=math.nan,
            total_observations=0,
            hyperparameters={"reference": dict(reference_spec or {})},
            ridge_lambda=ridge_lambda,
            final_training_period=FINAL_TRAINING_PERIOD_LABEL,
        )

    stats_by_year = {entry.year: entry for entry in stats}
    common_years = sorted(stats_by_year.keys())

    feature_dim = stats[0].S.shape[0]
    S_total = np.zeros((feature_dim, feature_dim), dtype=float)
    b_total = np.zeros(feature_dim, dtype=float)
    total_obs = 0

    for entry in stats:
        if entry.S.shape != (feature_dim, feature_dim):
            raise ValueError("Inconsistent feature dimension across yearly statistics.")
        S_total += entry.S
        b_total += entry.b
        total_obs += entry.n

    # Use vectorized accumulation
    total_sum_arr, total_count_arr, yearly_data_arr = _accumulate_climatology_arrays(stats)

    year_metrics: List[YearlyValidationMetrics] = []
    weighted_mse = 0.0
    weighted_mse_ref = 0.0
    accumulated_obs = 0

    for year in common_years:
        stat = stats_by_year[year]

        if total_obs - stat.n <= 0:
            continue

        S_excl = S_total - stat.S
        b_excl = b_total - stat.b
        beta = solve_normal_equations(S_excl, b_excl, ridge_lambda)

        y_pred = stat.X @ beta
        residuals_model = stat.y - y_pred

        # Vectorized reference calculation
        # Get excluded sums/counts
        year_sum_arr, year_count_arr = yearly_data_arr.get(year, (np.zeros_like(total_sum_arr), np.zeros_like(total_count_arr)))
        
        sum_excl_arr = total_sum_arr - year_sum_arr
        count_excl_arr = total_count_arr - year_count_arr
        
        # Map stat indices to array indices
        # Filter invalid indices first
        valid_indices = (
            (stat.utc_day_index >= 1)
            & (stat.utc_day_index <= 366)
            & (stat.utc_hour >= 0)
            & (stat.utc_hour < 24)
        )
        
        # Initialize ref_values with NaNs
        ref_values = np.full(stat.n, np.nan, dtype=float)
        
        # Only compute for valid indices
        if np.any(valid_indices):
            days = stat.utc_day_index[valid_indices]
            hours = stat.utc_hour[valid_indices]
            
            # Get counts for these points
            counts = count_excl_arr[days, hours]
            sums = sum_excl_arr[days, hours]
            
            # Where count > 0, compute mean
            valid_ref = counts > 0
            
            # We need to map back to the original `stat` array indices
            # valid_indices is a boolean mask of length stat.n
            # valid_ref is a boolean mask of length valid_indices.sum()
            
            # Create a mask for where we have a valid reference value
            # It must be a valid index AND have enough data in the baseline
            final_mask_subset = valid_ref
            
            # Calculate means where possible
            means = sums[valid_ref] / counts[valid_ref]
            
            # Assign back to ref_values
            # We need indices in the original array where valid_indices is True AND valid_ref is True
            # np.where(valid_indices)[0] gives indices where valid_indices is True
            # We take the subset of those where valid_ref is True
            target_indices = np.where(valid_indices)[0][valid_ref]
            ref_values[target_indices] = means

        finite_mask = np.isfinite(residuals_model) & np.isfinite(ref_values)
        if not np.any(finite_mask):
            continue

        model_errors = residuals_model[finite_mask]
        reference_errors = stat.y[finite_mask] - ref_values[finite_mask]

        mse_model = float(np.mean(np.square(model_errors)))
        mse_ref = float(np.mean(np.square(reference_errors)))
        rmse_model = math.sqrt(mse_model)
        skill = float("nan") if mse_ref <= 0.0 else float(1.0 - (mse_model / mse_ref))
        n_valid = int(model_errors.size)

        year_metrics.append(
            YearlyValidationMetrics(
                year=year,
                mse_model=mse_model,
                mse_ref=mse_ref,
                rmse=rmse_model,
                skill=skill,
                n=n_valid,
            )
        )
        weighted_mse += n_valid * mse_model
        weighted_mse_ref += n_valid * mse_ref
        accumulated_obs += n_valid

    if accumulated_obs == 0:
        global_rmse = math.nan
        global_skill = math.nan
    else:
        global_rmse = math.sqrt(weighted_mse / accumulated_obs)
        global_skill = float("nan") if weighted_mse_ref == 0.0 else float(
            1.0 - (weighted_mse / weighted_mse_ref)
        )

    hyperparameters = {
        "reference": dict(reference_spec or {}),
        "evaluation_time_base": "UTC",
        "model_time_base": "solar",
        "baseline": "climatology_mean per (utc_day, utc_hour), LOYO",
    }

    return LeaveOneYearOutReport(
        years=year_metrics,
        global_rmse=global_rmse,
        global_skill=global_skill,
        total_observations=accumulated_obs,
        hyperparameters=hyperparameters,
        ridge_lambda=ridge_lambda,
        final_training_period=FINAL_TRAINING_PERIOD_LABEL,
    )


__all__ = [
    "accumulate_climatology_maps",
    "evaluate_loyo",
]
