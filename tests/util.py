import numpy as np
import aero2tools

def assert_valid(actual, expected, margin = 0.02):
    rel_error = abs((actual - expected) / expected)
    assert (rel_error < margin)