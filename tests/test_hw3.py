from util import *

def test_problem2():
    
    T1 = Q_(300, 'K')
    P1 = Q_(1.2, 'atm')
    V1 = Q_(250, 'mps')
    
    state1 = Isen.from_velocity(T1, V1)
    state1.P = P1    
    
    actual1 = state1.T0.to('K').m
    actual2 = state1.Tstar.to('K').m
    actual3 = state1.P0.to('atm').m
    actual4 = state1.Pstar.to('atm').m
    actual5 = state1.mach_star.to('').m
    
    expected1 = 331.1
    expected2 = 275.9
    expected3 = 1.694
    expected4 = 0.895
    expected5 = 0.75
    
    assert_valid(actual1, expected1)
    assert_valid(actual2, expected2)
    assert_valid(actual3, expected3)
    assert_valid(actual4, expected4)
    assert_valid(actual5, expected5)

def test_problem3():
    
    P1 = Q_(1, 'atm')
    T1 = Q_(230, 'K')
    M1 = Q_(2.3, '')
    
    state1 = Isen(M1)
    state1.P = P1
    state1.T = T1
    
    actual1 = state1.T0.to('K').m
    actual2 = state1.P0.to('atm').m
    
    expected1 = 473.34
    expected2 = 12.5
    
    assert_valid(actual1, expected1)
    assert_valid(actual2, expected2)
    

def test_problem4():
    
    state1 = Isen(0.82)
    state2 = Isen(0.86)
    
    inter = IsenTranslate(state1, state2)
    
    inter.state1.P = Q_(1455.6, 'lbf/ft^2')
    inter.state1.T = Q_(483.04, 'rankine')
    
    actual1 = inter.state2.P.to('lbf/ft^2').m
    actual2 = inter.state2.T.to('rankine').m
    
    assert_valid(actual1, 1397)
    assert_valid(actual2, 477.0)
    