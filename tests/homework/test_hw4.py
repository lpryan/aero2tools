from tests.conftest import *

def test_problem1():
    Q_ = config.Q_
    
    P1 = Q_(1, 'atm')
    T1 = Q_(288, 'K')
    M1 = Q_(2.5, '')
    
    shock1 = Normal(M1, T = T1, P = P1)
    
    actual1 = shock1.P2.to('atm').m
    expected1 = 7.13
    assert_valid(actual1, expected1)
    
    actual2 = shock1.T2.to('K').m
    expected2 = 616
    assert_valid(actual2, expected2)
    
    actual3 = shock1.r2.to('kg/m^3').m
    expected3 = 4.086
    assert_valid(actual3, expected3)
    
    actual4 = shock1.mach2.to('').m
    expected4 = 0.513
    assert_valid(actual4, expected4)
    
    actual5 = shock1.P02.to('atm').m
    expected5 = 8.53
    assert_valid(actual5, expected5)
    
    actual6 = shock1.T02.to('K').m
    expected6 = 648
    assert_valid(actual6, expected6)
    
    actual7 = shock1.ds.to('J/(kg*K)').m
    expected7 = 199.5
    assert_valid(actual7, expected7)
    
def test_problem2():
    Q_ = config.Q_
    
    P1 = Q_(1, 'atm')
    P2 = Q_(10.33, 'atm')
    
    T2 = Q_(1350, 'rankine')
    
    shock1 = Normal()
    shock1.P2_P1 = P2 / P1
    
    shock1.P = P1
    shock1.P2 = P2
    shock1.T2 = T2
    
    actual = shock1.mach.to('').m
    expected = 3.0
    assert_valid(actual, expected)
    
    actual = shock1.T.to('rankine').m
    expected = 503.9
    assert_valid(actual, expected)
    
    actual = shock1.T02.to('rankine').m
    expected = 1411
    assert_valid(actual, expected)
    
    actual = shock1.P02.to('atm').m
    expected = 12.06
    assert_valid(actual, expected)
        
def test_problem3():
    Q_ = config.Q_
    
    P1 = Q_(0.1, 'atm')
    P02 = Q_(1.285, 'atm')
    
    shock1 = Normal()
    shock1.P02_P1 = P02 / P1
    
    actual = shock1.mach.to('').m
    expected = 3.1
    assert_valid(actual, expected)
    
def test_problem4():
    Q_ = config.Q_
    
    state1 = Isen(36, T = Q_(330, 'K'))
    
    actual = state1.T0.to('K').m
    expected = 85866
    assert_valid(actual, expected)

def test_problem5():
    Q_ = config.Q_
    
    shock1 = Normal(3, T = Q_(300, 'K'), P = Q_(1.5, 'atm'))
    
    actual = shock1.T0.to('K').m
    expected = 840
    assert_valid(actual, expected)
    
    actual = shock1.P0.to('atm').m
    expected = 55.10
    assert_valid(actual, expected)
    
    actual = shock1.r.to('kg/m^3').m
    expected = 1.765
    assert_valid(actual, expected)
    
    actual = shock1.mach2.to('').m
    expected = 0.475
    assert_valid(actual, expected)
    
    actual = shock1.T2.to('K').m
    expected = 803.70
    assert_valid(actual, expected)
    
    actual = shock1.T02.to('K').m
    expected = 840
    assert_valid(actual, expected)
    
    actual = shock1.P2.to('atm').m
    expected = 15.5
    assert_valid(actual, expected)
    
    actual = shock1.P02.to('atm').m
    expected = 18.091
    assert_valid(actual, expected)
    
    actual = shock1.r2.to('kg/m^3').m
    expected = 6.806
    assert_valid(actual, expected)
    
    actual = shock1.P02_P01.to('').m
    expected = 0.3283
    assert_valid(actual, expected)