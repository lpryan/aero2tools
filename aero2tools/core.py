import numpy as np 
from pint import UnitRegistry # type: ignore 


# --- Unit Setup --- 
ur = UnitRegistry() 
Q_ = ur.Quantity
ur.formatter.default_format = ".3f~P"

# --- Gas Constants --- 
R = Q_(287, "J/(kg*K)")


# --- GAMMA approximation ---
GAMMA = 1.4

CP = GAMMA * R / (GAMMA - 1)
CV = R / (GAMMA - 1)






# derivative approximation
def diff(func, x):
    epsilon = 1e-100   
    return np.imag(func(x + epsilon*1j)) / epsilon

def diff2(func, x):
    h = 1e-5
    
    f1 = diff(func, x)
    f2 = diff(func, x+h)
    
    return (f2 - f1)/h



# SPEED OF SOUND
class SpeedOfSound:
    @staticmethod
    def a(temp, GAMMA = GAMMA):
        return np.sqrt(GAMMA * R * temp)
    
    @staticmethod
    def mach(temp, velocity, GAMMA = GAMMA):
        return (velocity / np.sqrt(GAMMA * R * temp)).to('')
    
    @staticmethod
    def velocity(temp, mach, GAMMA = GAMMA):
        return SpeedOfSound.a(temp, GAMMA) * mach


# Ideal Gas Law
class IdealGas:
    @staticmethod
    def pressure(T, rho):
        return rho * R * T
    
    @staticmethod
    def temperature(P, rho):
        return P / (rho * R)
    
    @staticmethod
    def density(P, T):
        rho = P / (T * R)
        if isinstance(rho, Q_): rho.ito("kg/m^3")
        return rho