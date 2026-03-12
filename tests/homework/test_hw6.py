from tests.conftest import *

def test_problem1a(subtests):
    Q_ = config.Q_
    
    M1 = Q_(2.2, '')
    theta1 = Q_(12.5, 'deg')
    
    shock1 = Oblique.IsenTheta(M1, theta1)
    
    assert_subtest(subtests, 'i', shock1.theta_max.to('deg').m, 26.10)
    assert_subtest(subtests, 'ii', shock1.mu.to('deg').m, 27.04)    
    assert_subtest(subtests, 'iii', shock1.beta.to('deg').m, 38.441)

def test_problem1b():
    Q_ = config.Q_
    
    M1 = Q_(2, '')
    theta1 = Q_(25, 'deg')
    
    try:
        shock1 = Oblique.IsenTheta(M1, theta1)
    
    except ValueError:
        assert True
        return
    
    assert False
    

def test_problem2(subtests):
    Q_ = config.Q_
    
    shock1 = Oblique.IsenTheta(2.5, Q_(22.5, 'deg'), P=Q_(2, 'atm'), T=Q_(260, 'K'))
    
    assert_subtest(subtests, 'a', shock1.beta.to('deg').m, 46)
    assert_subtest(subtests, 'b', shock1.P2.to('atm').m, 7.29)
    assert_subtest(subtests, 'c', shock1.T2.to('K').m, 399.78)
    assert_subtest(subtests, 'd', shock1.mach2.to('').m, 1.52)


def test_problem3():
    Q_ = config.Q_
    
    state0 = Isen(3.5, P = Q_(0.5, 'atm'))
    
    shocks = IsenTracker(state0)
    
    shocks.addShock(Oblique, theta = Q_(30.2, 'deg'))
    
    shocks.P2
    
    shocks.addShock(Normal)
    
    assert_valid(shocks.P03.to('atm').m, 15.39)
    
    
def test_problem4():
    Q_ = config.Q_
    
    state0 = Isen(4, P = Q_(2, 'atm'))
    shocks = IsenTracker(state0)
    
    shocks.addShock(Normal)
    
    dP = shocks.P01 - shocks.P02
    
    assert_valid(dP.to('atm').m, 261.5)

def test_problem5():
    Q_ = config.Q_
    
    state0 = Isen(4, P = Q_(2, 'atm'))
    shocks = IsenTracker(state0)
    
    shocks.addShock(Oblique, theta = Q_(25.3, 'deg'))
    shocks.addShock(Normal)
    
    
    dP = shocks.P01 - shocks.P03
    
    assert_valid(dP.to('atm').m, 208.3)
    
def test_problem6():
    Q_ = config.Q_
    
    state0 = Isen(4, P = Q_(2, 'atm'))
    shocks = IsenTracker(state0)
    
    shocks.addShock(Oblique, theta = Q_(25.3, 'deg'))
    shocks.addShock(Oblique, theta = Q_(20, 'deg'))
    shocks.addShock(Normal)
    
    dP = shocks.P01 - shocks.P04
    
    assert_valid(dP.to('atm').m, 176.2)   


def test_problem7(subtests):
    Q_ = config.Q_
    
    state0 = Isen(3.2, P = Q_(2, 'atm'), T = Q_(520, 'rankine'))
    shocks = IsenTracker(state0)
    
    shocks.addShock(Oblique, theta = Q_(18.2, 'deg'))
    shocks.addShock(Oblique, theta = Q_(18.2, 'deg'))
        
    assert_subtest(subtests, 'a', shocks.mach3.to('').m, 1.506)
    assert_subtest(subtests, 'b', (shocks.beta2 - shocks.theta2).to('deg').m, 26.75)
    assert_subtest(subtests, 'c', shocks.P3.to('atm').m, 19.63)
    assert_subtest(subtests, 'd', shocks.T3.to('rankine').m, 1090.36)    