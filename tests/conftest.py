from aero2tools import *
import pytest

def assert_valid(actual, expected, margin = 0.02):
    assert actual == pytest.approx(expected, rel = margin)
    
def assert_subtest(subtests, num, actual, expected, margin = 2e-2):
    
    with subtests.test(num):
        assert_valid(actual, expected, margin)