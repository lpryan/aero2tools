import numpy as np
import functools
import inspect
import uuid

from pint import UnitRegistry
import pint

# ===========================
# Default Config
# ===========================

class Config:
    
    def __init__(self):
        
        self.ur = UnitRegistry()
        self.Q_ = self.ur.Quantity
        self.ur.formatter.default_format = ">6.4g~P"
        
        self.R = self.Q_(287, "J/(kg*K)")
        self.GAMMA = 1.4
        
    @property
    def CP(self):
        return self.GAMMA * self.R / (self.GAMMA - 1)
    
    @property
    def CV(self):
        return self.R / (self.GAMMA - 1)
    
    
    def wrap(self, out_unit, in_units, strict = True):
        pint_decorator = self.ur.wraps(out_unit, in_units, strict)
        def decorator(func):
            wrapped = pint_decorator(func)
            functools.update_wrapper(wrapped, func)
            wrapped.__signature__ = inspect.signature(func)
            return wrapped
        return decorator

config = Config()



# ===========================
# Ideal Gas Law
# ===========================

class IdealGas:
    """
    Ideal‑gas relations for computing missing state variables.

    Provides static methods to compute any one of:
        - (pres) Pressure  [atm]
        - (temp) Temperature [K]
        - (dens) Density [kg/m^3]

    given the other two, using the ideal gas law:
        P = ρ R T
    """
    
    @staticmethod
    @config.wrap('atm', ('K', 'kg/m^3'), False)
    def pres(temp: pint.Quantity, dens: pint.Quantity) -> pint.Quantity:
        """
        Pressure from temperature and density.\n
        temp [K], dens [kg/m^3] -> pressure [atm]
        """
        pres = dens * config.R.m * temp / 101325
        return pres
    
    @staticmethod
    @config.wrap('K', ('atm', 'kg/m^3'), False)
    def temp(pres: pint.Quantity, dens: pint.Quantity) -> pint.Quantity:
        """
        temperature from pressure and density.\n
        pres [atm], dens [kg/m^3] -> temperature [K]
        """
        temp = 101325 * pres / (dens * config.R.m)
        return temp
    
    @staticmethod
    @config.wrap('kg/m^3', ('atm', 'K'), False)
    def dens(pres: pint.Quantity, temp: pint.Quantity) -> pint.Quantity:
        """
        Density from pressure and temperature.\n
        pres [atm], temp [K] -> density [kg/m^3]
        """
        dens = 101325 * pres / (temp * config.R.m)
        return dens


# ===========================
# Speed of Sound
# ===========================

class SpeedOfSound:
    """
    Relations for computing aspects of the speed of sound.
    
    Provides static methods to compute any one of:
        - (a) Speed of sound [m/s]
        - (mach) Mach []
        - (vel) Velocity [m/s]
        
    provided temperature and velocity or mach
    """
    
    @staticmethod
    @config.wrap('m/s', ('K', None), False)
    def a(temp, GAMMA = config):
        """
        Speed of sound from temperature.\n
        temp [K] -> speed of sound [m/s]
        """
        if GAMMA is config: GAMMA = config.GAMMA
        return np.sqrt(GAMMA * config.R.m * temp)
    
    @staticmethod
    @config.wrap('', ('K', 'm/s', None), False)
    def mach(temp, vel, GAMMA = config):
        """
        Mach from temperature and velocity.\n
        temp [K], vel [m/s] -> mach []
        """
        if GAMMA is config: GAMMA = config.GAMMA
        return vel / SpeedOfSound.a(temp, GAMMA).m
    
    @staticmethod
    @config.wrap('m/s', ('K', '', None), False)
    def vel(temp, mach, GAMMA = config):
        """
        Speed of sound from temperature.\n
        temp [K], mach [] -> velocity [m/s]
        """
        if GAMMA is config: GAMMA = config.GAMMA
        return SpeedOfSound.a(temp, GAMMA).m * mach
