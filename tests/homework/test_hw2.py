from tests.conftest import *

def test_problem2():
    Q_ = config.Q_    
    
    Pinf = Q_(0.61, 'atm')
    rinf = Q_(0.819, 'kg/m^3')
    Vinf = Q_(300, 'mps')
    
    P2 = Q_(0.430, 'atm')
    
    Tinf = IdealGas.temp(Pinf, rinf)
    state1 = Isen.from_vel(Tinf, Vinf)
    state1.P = Pinf
    
    inter = IsenTranslate.from_pres1(state1, P2 = P2)
    
    V2 = SpeedOfSound.vel(inter.T2, inter.mach2)
    
    actual = V2.to('mps').m
    expected = 374.470
    
    assert_valid(actual, expected)

def test_problem3():
    Q_ = config.Q_
    
    Pinf = Q_(0.61, 'atm')
    rinf = Q_(0.819, 'kg/m^3')
    Vinf = Q_(300, 'mps')
    
    P2 = Q_(0.540, 'atm')

    Tinf = IdealGas.temp(Pinf, rinf)
    state1 = Isen.from_vel(Tinf, Vinf)
    
    inter = IsenTranslate.from_pres1(state1, P1 = Pinf, P2 = P2)
    
    V2 = SpeedOfSound.vel(inter.T2, inter.mach2)
    
    actual = V2.to('mps').m
    expected = 328.757
    
    assert_valid(actual, expected)

def test_problem4():
    Q_ = config.Q_
    
    Pinf = Q_(0.61, 'atm')
    rinf = Q_(0.819, 'kg/m^3')
    Vinf = Q_(300, 'mps')
    
    P2 = Q_(0.3, 'atm')
    
    Tinf = IdealGas.temp(Pinf, rinf)
    state1 = Isen.from_vel(Tinf, Vinf)
    
    inter = IsenTranslate.from_pres1(state1, P1 = Pinf, P2 = P2)
    
    V2 = inter.vel2
    
    actual = V2.to('mps').m
    expected = 432.382
    
    assert_valid(actual, expected)
    
    
    
    
    
    
    
    
    
    
    
    
    
    




