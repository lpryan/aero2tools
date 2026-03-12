from tests.conftest import *

def test_betamax():
    Q_ = config.Q_
    
    shock1 = Oblique(2.3, Q_(1, 'rad'))
    
    actual = shock1.beta_max.to('rad').m
    expected = 1.12841
    
    assert_valid(actual, expected)
    
def test_thetamax():
    Q_ = config.Q_
    
    shock1 = Oblique(3, Q_(50, 'deg'))
    
    actual = shock1.theta_max.to('deg').m
    expected = 34.07344
    
    assert_valid(actual, expected)


def test_mu1():
    Q_ = config.Q_
    
    shock1 = Isen(3)
    
    actual = shock1.mu.to('deg').m
    expected = 19.47
    
    assert_valid(actual, expected)
    
def test_mu2():
    Q_ = config.Q_
    
    shock1 = Normal(3)
    
    actual = shock1.mu.to('deg').m
    expected = 19.47
    
    assert_valid(actual, expected)

def test_mu3():
    Q_ = config.Q_
    
    shock1 = Oblique(3, Q_(50, 'deg'))
    
    actual = shock1.mu.to('deg').m
    expected = 19.47
    
    assert_valid(actual, expected)
    

def test_nu1():
    Q_ = config.Q_
    
    state0 = Isen(2.5)
    assert_valid(state0.nu.to('deg').m, 39.1236)

def test_nu2():
    Q_ = config.Q_
    
    state0 = Isen(1.35)
    assert_valid(state0.nu.to('deg').m, 7.561)



def test_propagation():
    
    Q_ = config.Q_
    
    state0 = Isen(3.5, P = Q_(0.5, 'atm'))
    
    shock1 = Oblique.IsenBeta(state0, Q_(60, 'deg'))
    assert_valid(state0.P0.to('atm').m, shock1.P0.to('atm').m)
    
    shock2 = Oblique.IsenBeta(shock1, Q_(60, 'deg'))
    assert_valid(shock1.P02.to('atm').m, shock2.P0.to('atm').m)
    
    shock3 = Normal(shock2)
    assert_valid(shock2.P02.to('atm').m, shock3.P0.to('atm').m)
    

def test_tracker():
    
    Q_ = config.Q_
    
    state0 = Isen(3.5, P = Q_(0.5, 'atm'))
    
    # manual propagation
    shock1 = Oblique.IsenTheta(state0, Q_(30, 'deg'))
    assert_valid(state0.P0.to('atm').m, shock1.P0.to('atm').m)
    
    shock2 = Oblique.IsenBeta(shock1, Q_(60, 'deg'))
    assert_valid(shock1.P02.to('atm').m, shock2.P0.to('atm').m)
    
    shock3 = Normal(shock2)
    assert_valid(shock2.P02.to('atm').m, shock3.P0.to('atm').m)
    
    # automatic propagation
    shocks = IsenTracker(state0)
    assert_valid(shocks.P.to('atm').m, state0.P.to('atm').m)
    
    shocks.addShock(Oblique, theta = Q_(30, 'deg'))
    assert_valid(shocks.P2.to('atm').m, shock1.P2.to('atm').m)
    
    shocks.addShock(Oblique, beta = Q_(60, 'deg'))
    assert_valid(shocks.P3.to('atm').m, shock2.P2.to('atm').m)
    
    shocks.addShock(Normal)
    assert_valid(shocks.P4.to('atm').m, shock3.P2.to('atm').m)

