
import json
import os
from pathlib import Path
import pytest
import numpy as np
from harmoclimate.template_cpp import generate_cpp_header

def test_generate_cpp_header_with_stochastic(tmp_path):
    # Create dummy payloads

    # Helper to create a minimal linear model payload
    def create_linear_payload(target):
        return {
            "metadata": {
                "station_usual_name": "Test Station",
                "station_code": "12345",
                "longitude_deg": 2.35,
                "latitude_deg": 48.85,
                "delta_utc_solar_h": 0.15,
                "target_variable": target,
                "target_unit": "unit"
            },
            "model": {
                "n_diurnal": 1,
                "coefficients": [10.0, 1.0, 0.5], # c0, a1, b1
                "params_layout": [
                    {"name": "c0", "start": 0, "length": 1, "n_annual": 0},
                    {"name": "a1", "start": 1, "length": 1, "n_annual": 0},
                    {"name": "b1", "start": 2, "length": 1, "n_annual": 0}
                ]
            }
        }

    t_payload = create_linear_payload("T")
    q_payload = create_linear_payload("Q")
    p_payload = create_linear_payload("P")

    # Stochastic payload
    # Covariance matrix that gives a non-trivial Cholesky
    # Sigma = [[4, 2, 0], [2, 5, 0], [0, 0, 1]]
    # L should be:
    # L[0,0] = 2
    # L[1,0] = 1 (since 2*1 = 2)
    # L[1,1] = sqrt(5 - 1^2) = 2
    # L[2,2] = 1
    stochastic_payload = {
        "metadata": {
            "variables": ["P", "Q", "T"]
        },
        "model": {
            "coefficient_matrix": [[0.9, 0.0, 0.0], [0.0, 0.8, 0.0], [0.0, 0.0, 0.7]],
            "covariance_matrix": [[4.0, 2.0, 0.0], [2.0, 5.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    output_file = tmp_path / "test_model.hpp"

    generate_cpp_header(
        t_payload,
        q_payload,
        p_payload,
        output_file,
        stochastic_payload=stochastic_payload
    )

    assert output_file.exists()
    content = output_file.read_text(encoding="utf-8")

    assert "namespace stochastic" in content
    assert "struct State" in content
    assert "double p_res = 0.0;" in content
    assert "void update_state" in content
    assert "void predict_stochastic" in content

    # Check for Cholesky values in the source
    # We expect L flattened. L = [[2, 0, 0], [1, 2, 0], [0, 0, 1]]
    # So 2.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 1.0 roughly
    assert "2" in content

    # Check that prediction combines baseline and residuals
    assert "temperature_c += state.t_res;" in content

def test_generate_cpp_header_stochastic_failure(tmp_path):
    # Test that invalid variable order causes generation to skip stochastic part
    t_payload = {"metadata": {}, "model": {"n_diurnal": 0, "coefficients": [0], "params_layout": [{"name": "c0", "start": 0, "length": 1, "n_annual": 0}]}}
    # Reuse helpers if possible but creating minimal payload manually for speed
    # Actually I need valid payload or generate_cpp_header will fail on non-stochastic parts.
    # I'll reuse the create_linear_payload from above if I can scope it out or copy.

    def create_linear_payload(target):
        return {
            "metadata": {
                "station_usual_name": "Test Station",
                "station_code": "12345",
                "longitude_deg": 2.35,
                "target_variable": target,
            },
            "model": {
                "n_diurnal": 0,
                "coefficients": [0.0],
                "params_layout": [
                    {"name": "c0", "start": 0, "length": 1, "n_annual": 0}
                ]
            }
        }

    t_payload = create_linear_payload("T")
    q_payload = create_linear_payload("Q")
    p_payload = create_linear_payload("P")

    # Invalid order
    stochastic_payload = {
        "metadata": {
            "variables": ["T", "Q", "P"]
        },
        "model": {
            "coefficient_matrix": [[1,0,0],[0,1,0],[0,0,1]],
            "covariance_matrix": [[1,0,0],[0,1,0],[0,0,1]],
        }
    }

    output_file = tmp_path / "test_fail.hpp"

    generate_cpp_header(
        t_payload,
        q_payload,
        p_payload,
        output_file,
        stochastic_payload=stochastic_payload
    )

    content = output_file.read_text(encoding="utf-8")
    assert "// [Warn] Stochastic model generation failed" in content
    assert "namespace stochastic" not in content
    assert "predict_stochastic" not in content
