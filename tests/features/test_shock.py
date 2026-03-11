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