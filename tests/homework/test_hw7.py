from tests.conftest import *

def test_problem2(subtests):
    Q_ = config.Q_
    
    state0 = Isen(2, P = Q_(0.6, 'atm'), T = Q_(590, 'rankine'))
    shocks = IsenTracker(state0)

    shocks.addShock(Expansion, theta = Q_(23.38, 'deg'))
    
    assert_subtest(subtests, 'a', shocks.mach2.to('').m, 3)
    
    assert_subtest(subtests, 'b', shocks.P2.to('atm').m, 0.128)
    
    assert_subtest(subtests, 'c', shocks.T2.to('rankine').m, 379.265)
    
    assert_subtest(subtests, 'd', shocks.r2.to('slug/ft^3').m, 4.154e-4)
    
    assert_subtest(subtests, 'e', shocks.P02.to('atm').m, 4.695)
    
    assert_subtest(subtests, 'f', shocks.T02.to('rankine').m, 1062)
    
    assert_subtest(subtests, 'g', shocks.fwd1.to('deg').m, 30)

    assert_subtest(subtests, 'h', shocks.rwd1.to('deg').m, -3.91)
    


def test_problem3():
    Q_ = config.Q_
    
    state0 = Isen(1.58, P = Q_(1.25, 'atm'))
    shocks = IsenTracker(state0)
    
    shocks.addShock(Expansion, P2 = Q_(0.1306, 'atm'))
    
    assert_valid(shocks.theta1.to('deg').m, 36.415)


def test_problem4(subtests):
    Q_ = config.Q_
    
    state0 = Isen(3, T = Q_(260, 'K'), P = Q_(1, 'atm'))
    shocks = IsenTracker(state0)
    
    shocks.addShock(Oblique, theta = Q_(30.6, 'deg'))
    shocks.addShock(Expansion, theta = Q_(30.6, 'deg'))
    
    
    assert_subtest(subtests, 'a', shocks.mach3.to('').m, 2.475)
    assert_subtest(subtests, 'b', shocks.P3.to('atm').m, 1.206)
    assert_subtest(subtests, 'c', shocks.T3.to('K').m, 327.162)
    
    
    