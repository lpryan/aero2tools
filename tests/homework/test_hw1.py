from tests.conftest import *




def test_problem1():
    Q_ = config.Q_    
    
    T = Q_(968, 'rankine')
    P = Q_(7.8, 'atm')
    
    actual = IdealGas.dens(P, T).to('slug/ft^3').m
    expected = 0.0099
    
    assert_valid(actual, expected)


def test_problem2():
    Q_ = config.Q_
    
    T = Q_(936, 'rankine')
    P = Q_(7.8, 'atm')
    
    actual1 = config.CP.to('ft*lbf/slug/rankine').m
    actual2 = config.CV.to('ft*lbf/slug/rankine').m
    actual3 = (config.CV * T).to('ft*lbf/slug').m
    actual4 = (config.CP * T).to('ft*lbf/slug').m
    
    expected1 = 6008
    expected2 = 4291
    expected3 = 4.017e6
    expected4 = 5.623e6
    
    assert_valid(actual1, expected1)
    assert_valid(actual2, expected2)
    assert_valid(actual3, expected3)
    assert_valid(actual4, expected4)


def test_problem3():
    Q_ = config.Q_
    
    T = Q_(519, 'rankine')
    P = Q_(7.8, 'atm')
    
    actual1 = config.CP.to('ft*lbf/slug/rankine').m
    actual2 = config.CV.to('ft*lbf/slug/rankine').m
    actual3 = (config.CV * T).to('ft*lbf/slug').m
    actual4 = (config.CP * T).to('ft*lbf/slug').m
    
    expected1 = 6008
    expected2 = 4291
    expected3 = 2.227e6
    expected4 = 3.118e6
    
    assert_valid(actual1, expected1)
    assert_valid(actual2, expected2)
    assert_valid(actual3, expected3)
    assert_valid(actual4, expected4)