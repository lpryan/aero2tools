from util import *

def test_problem2():
    
    Pinf = Q_(0.61, 'atm')
    rinf = Q_(0.819, 'kg/m^3')
    Vinf = Q_(300, 'mps')
    
    P2 = Q_(0.430, 'atm')
    
    Tinf = IdealGas.temperature(Pinf, rinf)
    state1 = Isen(SpeedOfSound.mach(Tinf, Vinf))
    
    inter = IsenTranslate.from_pressure1(state1, P1 = Pinf, P2 = P2)
    inter.state1.T = Tinf
    inter.state1.P = Pinf
    
    T2 = inter.state2.T
    
    V2 = SpeedOfSound.velocity(T2, inter.state2.mach)
    
    actual = V2.to('mps').m
    expected = 374.470
    
    assert_valid(actual, expected)

def test_problem3():
    
    Pinf = Q_(0.61, 'atm')
    rinf = Q_(0.819, 'kg/m^3')
    Vinf = Q_(300, 'mps')
    
    P2 = Q_(0.540, 'atm')
    
    Tinf = IdealGas.temperature(Pinf, rinf)
    mach1 = SpeedOfSound.mach(Tinf, Vinf)
    
    state1 = Isen(mach1)
    inter = IsenTranslate.from_pressure1(state1, P1 = Pinf, P2 = P2)
    
    inter.state1.P = Pinf
    inter.state1.T = Tinf
    
    mach2 = inter.state2.mach
    T2 = inter.state2.T
    
    V2 = SpeedOfSound.velocity(T2, mach2)
    
    
    actual = V2.to('mps').m
    expected = 328.757
    
    assert_valid(actual, expected)

def test_problem4():
    
    Pinf = Q_(0.61, 'atm')
    rinf = Q_(0.819, 'kg/m^3')
    Vinf = Q_(300, 'mps')
    
    P2 = Q_(0.3, 'atm')
    
    Tinf = IdealGas.temperature(Pinf, rinf)
    state1 = Isen.from_velocity(Tinf, Vinf)
    
    inter = IsenTranslate.from_pressure1(state1, P1 = Pinf, P2 = P2)
    
    inter.state1.P = Pinf
    inter.state1.T = Tinf
        
    V2 = inter.state2.vel
    
    actual = V2.to('mps').m
    expected = 432.382
    
    assert_valid(actual, expected)
