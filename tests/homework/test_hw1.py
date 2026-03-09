from tests.conftest import *




def test_problem1():
    Q_ = config.Q_    
    
    T = Q_(968, 'rankine')
    P = Q_(7.8, 'atm')
    
    actual = IdealGas.dens(P, T).to('slug/ft^3').m
    expected = 0.0099
    
    assert_valid(actual, expected)