from tests.conftest import *


def test_problem2(subtests):
    
    Q_ = config.Q_
    
    shock1 = Oblique(4, Q_(30, 'deg'), P = Q_(2.65e4, 'N/m^2'), T = Q_(223.3, 'K'))
    
    assert_subtest(subtests, "a", shock1.P2.to('N/m^2').m, 1.19e5)
    
    assert_subtest(subtests, "b", shock1.T2.to('K').m, 376.82)
    
    assert_subtest(subtests, "c", shock1.mach2.to('').m, 2.73)
    
    assert_subtest(subtests, "d", shock1.P02.to('N/m^2').m, 2.9e6)
    
    assert_subtest(subtests, "e", shock1.T02.to('K').m, 937.86)
    
    assert_subtest(subtests, "f", shock1.ds.to('J/(kg*K)').m, 93.93)
    

def test_problem4():
    Q_ = config.Q_
    
    shock1 = Oblique(3, Q_(36.87, 'deg'), P = Q_(1, 'atm'))
    
    assert shock1.P02.to('atm').m == pytest.approx(29.85, rel=2e-2)

def test_problem7():
    Q_ = config.Q_
    
    shock1 = Oblique(2.4, Q_(34, 'deg'))
    
    assert shock1.theta.to('deg').m == pytest.approx(11, rel=2e-2)