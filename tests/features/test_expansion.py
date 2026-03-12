from tests.conftest import *

def test_expansion1(subtests):
    Q_ = config.Q_
    
    state0 = Isen(2)
    exp1 = Expansion.from_dtheta1(state0, Q_(12, 'deg'))
    
    assert_subtest(subtests, 'a', exp1.mach2.to('').m, 2.468)
    
    assert_subtest(subtests, 'b', exp1.P2_P1.to('').m, 0.4811)
    
    assert_subtest(subtests, 'c', exp1.T2_T1.to('').m, 0.8114)
    
def test_expansion2(subtests):
    Q_ = config.Q_
    
    state0 = Isen(1.1, T = Q_(273, 'K'), P = Q_(1.2, 'atm'))
    exp1 = Expansion.from_dtheta1(state0, Q_(40, 'deg'))
    
    assert_subtest(subtests, 'a', exp1.mach2.to('').m, 2.6)
    assert_subtest(subtests, 'b', exp1.P2.to('atm').m, 0.129)
    
    assert_subtest(subtests, 'c', exp1.r1.to('kg/m^3').m, 1.551)
    assert_subtest(subtests, 'd', exp1.r2.to('kg/m^3').m, 0.316)

def test_expansion_tracker(subtests):
    Q_ = config.Q_
    
    state0 = Isen(1.1, T = Q_(273, 'K'), P = Q_(1.2, 'atm'))
    shocks = IsenTracker(state0)
    
    shocks.addShock(Expansion, theta = Q_(40, 'deg'))
    
    assert_subtest(subtests, 'a', shocks.mach2.to('').m, 2.6)
    assert_subtest(subtests, 'b', shocks.P2.to('atm').m, 0.129)
    
    assert_subtest(subtests, 'c', shocks.r1.to('kg/m^3').m, 1.551)
    assert_subtest(subtests, 'd', shocks.r2.to('kg/m^3').m, 0.316)

def test_expansion_tracker2(subtests):
    Q_ = config.Q_
    
    state0 = Isen(1.1, T = Q_(273, 'K'), P = Q_(1.2, 'atm'))
    shocks = IsenTracker(state0)
    
    shocks.addShock(Expansion, P2 = Q_(0.129, 'atm'))
    
    assert_subtest(subtests, 'a', shocks.mach2.to('').m, 2.6)
    assert_subtest(subtests, 'b', shocks.T2.to('K').m, 144.383)
    
    assert_subtest(subtests, 'c', shocks.dtheta1.to('deg').m, 40)