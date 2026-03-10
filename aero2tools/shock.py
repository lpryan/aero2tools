from .core import *
from .isentropic import *

import numpy as np
import re


# ============================================================ 
# General Shock
# ============================================================

class Shock(Isen):
    
    def __init__(self, M, **kwargs):
        if (M < 1): raise ValueError(f"shock cannot occure with subsonic mach (M = {M})")
        
        if isinstance(M, Normal) or isinstance(M, Oblique):
            M = M.state2
            
        if isinstance(M, Isen):
            state = M
            mach = config.Q_(state.mach).m
            
            super().__init__(mach, **kwargs)
            self.propagate_state(state)
            
        else:
            super().__init__(M, **kwargs)
    
    def propagate_state(self, state: Isen):
        return super().propagate_state(state)
    
    def propagate_state2(self, state: 'Shock'):
        
        self.mach = state.mach2
    
        try: self.T = state.T2
        except TypeError: pass
        
        try: self.P = state.P2
        except TypeError: pass
        
        try: self.r = state.r2
        except TypeError: pass
    
    def __getattr__(self, name):
        
        try:
            match name:
                case 'T2': 
                    return self.T2_T1 * self.T
                
                case 'P2': 
                    return self.P2_P1 * self.P
                
                case 'r2': 
                    return self.r2_r1 * self.r
                
                case 'P02': 
                    return self.P02_P01 * self.P0
                
                case 'T02': 
                    return self.T0
                
                case 'ds':
                    return -config.R * np.log(self.P02_P01)
                
                case 'T0_T2':
                    return self.T0_T / self.T2_T1
                
                case 'P0_P2':
                    return self.P0_P / self.P2_P1
                
                case 'r0_r2':
                    return self.r0_r / self.r2_r1
        
        except TypeError:
            raise TypeError(f"Invalid Type: please assign values to normal shock")
        
        raise AttributeError(f"Normal has no attirubte \'{name}\'")
    

# ============================================================ 
# Normal Shock
# ============================================================

class Normal(Shock):
    
    def __init__(self, M = 1, **kwargs):
        super().__init__(M, **kwargs)
        
    
    @property
    def state2(self):
        return self._state2
    
    # --------------
    @property
    def mach2(self):
        
        M2_num = np.pow(super().mach, 2) * (self.GAMMA - 1) + 2
        M2_den = 2*self.GAMMA*np.pow(super().mach, 2) - (self.GAMMA - 1)
        
        return np.sqrt(M2_num / M2_den)
    
    @mach2.setter
    def mach2(self, mach2):
        
        num = np.sqrt((self.GAMMA - 1)*np.pow(mach2, 2) + 2)
        den = np.sqrt(self.GAMMA * (2 * np.pow(mach2, 2) - 1) + 1)
        
        self.mach = num / den
        
        state = Shock(self.mach2)
        state.propagate_state2(self)
        self._state2 = state    
    
    # --------------
    @property
    def r2_r1(self):
        
        num = (self.GAMMA + 1) * np.pow(super().mach, 2)
        den = (self.GAMMA - 1)*np.pow(super().mach, 2) + 2
        
        return num / den
    
    @r2_r1.setter
    def r2_r1(self, r2_r1):
        
        num = 2 * r2_r1
        den = (self.GAMMA + 1) - r2_r1 * (self.GAMMA - 1)
        
        M = np.sqrt(num / den)
        self.mach = M
    
    # --------------
    @property
    def P2_P1(self):
        
        num = 2*self.GAMMA*np.pow(super().mach, 2) - (self.GAMMA - 1)
        den = self.GAMMA + 1
        
        return num / den
    
    @P2_P1.setter
    def P2_P1(self, P2_P1):
        
        num = P2_P1 * (self.GAMMA + 1) + (self.GAMMA - 1)
        den = 2 * self.GAMMA
        
        M = np.sqrt(num / den)
        self.mach = M
    
    # --------------
    @property
    def T2_T1(self):
        return self.P2_P1 / self.r2_r1
    
    @T2_T1.setter
    def T2_T1(self, T2_T1):
        if (T2_T1 <= 1): raise ValueError("T2/T1 must be greater than 1 for Normal Shocks")
        
        g = self.GAMMA
        x = T2_T1
        
        a = ( g**2 * x**2 + 2*g**2 * x + g**2 + 2*g * x**2 - 12*g * x + 2*g + x**2 + 2*x + 1 ) 
        
        # Full numerator and denominator of the Mach expression 
        num = (g + 1) * np.sqrt(a) + (g**2 + 1) * (x + 1) + 2*g*x - 6*g 
        den = 4 * (g - 1) * g 
        
        M = np.sqrt(num / den) 
        self.mach = M
    
    # --------------
    @property
    def P02_P01(self):
        return self.P2_P1 * np.pow(self.T2_T1, -(self.GAMMA)/(self.GAMMA - 1))
    
    @P02_P01.setter
    def P02_P01(self, P02_P01):
        pass
    
    # --------------
    @property
    def P02_P1(self):
        M = super().mach
        
        t1 = np.pow((self.GAMMA + 1)*M, 2)/(4*self.GAMMA*np.pow(M, 2) - 2*(self.GAMMA - 1))
        t2 = ((1 - self.GAMMA) + 2*self.GAMMA*np.pow(M,2)) / (self.GAMMA + 1)
        
        return np.pow(t1, self.GAMMA/(self.GAMMA - 1)) * t2
    
    @P02_P1.setter
    def P02_P1(self, P02_P1):
        
        num = 9*np.pow(self.GAMMA + 1, 2)
        den = 2*(17*self.GAMMA + 1)        
        
        K = np.pow(num/den, self.GAMMA / (self.GAMMA - 1))
        
        
        a_num = 2*self.GAMMA*(289*self.GAMMA + 19)
        a_den = 9*(self.GAMMA + 1)*(17*self.GAMMA+1)
        a1 = a_num / a_den * K
        
        b1 = (34*self.GAMMA)/(3*(self.GAMMA+1)) * K
        
        c1 = (17*self.GAMMA + 1)/(self.GAMMA + 1) * K
        
        a = a1 / 2
        b = b1 - 3*a1
        c = 9/2*a1 - 3*b1 + c1
        
        M = (-b + np.sqrt(np.pow(b, 2) - 4*a*(c - P02_P1)))/2/a
        self.mach = M

    # --------------
    def __setattr__(self, name, value):
        
        match name:
            case 'T2':
                self.T = value / self.T2_T1
                return
            
            case 'P2':
                self.P = value / self.P2_P1
                return
                
            case 'r2':
                self.r = value / self.r2_r1
                return
        
        return super().__setattr__(name, value)
        

# ============================================================ 
# Oblique Shock
# ============================================================

def thetaMachBeta(mach1, beta):
    theta = beta - np.atan((5/(np.pow(mach1*np.sin(beta), 2)) + 1) * np.tan(beta)/6)
    return theta

class Oblique(Shock):
    
    def __init__(self, M, beta, **kwargs):
        super().__init__(M, **kwargs)
        
        #self.beta = beta
        
    @property
    def beta_max(self):
        mach = self.mach.to('').m
        beta_max = optimize.optimize(lambda B: thetaMachBeta(mach, B), np.pi / 4)
        return config.Q_(beta_max,'rad')
        
    @property
    def theta_max(self):
        return thetaMachBeta(self.mach, self.beta_max)
    
    
    













