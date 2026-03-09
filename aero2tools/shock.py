from .core import *
from .isentropic import *

import numpy as np
import re

class Shock(Isen):
    
    def __init__(self, M, **kwargs):
        if (M < 1): raise ValueError(f"shock cannot occur with subsonic mach (M = {M})")
        super().__init__(M, **kwargs)
    
    
    def propagate_state(self, state: Isen):
    
        self.mach = state.mach
    
        try: self.T = state.T
        except TypeError: pass
        
        try: self.P = state.P
        except TypeError: pass
        
        try: self.rho = state.rho
        except TypeError: pass
    
    
    def propagate_state2(self, state: 'Shock'):
        
        self.mach = state.mach2
    
        try: self.T = state.T2
        except TypeError: pass
        
        try: self.P = state.P2
        except TypeError: pass
        
        try: self.rho = state.rho2
        except TypeError: pass


# ============================================================ 
# Normal Shock
# ============================================================

class Normal(Shock):
    
    def __init__(self, M = 1, **kwargs):
                
        if isinstance(M, Normal) or isinstance(M, Oblique):
            M = M.state2
        
        if isinstance(M, Isen):
            state = M
            mach = Q_(state.mach).m
        
            super().__init__(mach, **kwargs)
            self.propagate_state(state)
        
        else:
            super().__init__(M, **kwargs)
            
        
            
            
    
    @property
    def state2(self):
        return self._state2
    
    @property
    def mach2(self):
        
        M2_num = np.pow(super().mach, 2) * (GAMMA - 1) + 2
        M2_den = 2*GAMMA*np.pow(super().mach, 2) - (GAMMA - 1)
        
        return np.sqrt(M2_num / M2_den)
    
    @mach2.setter
    def mach2(self, M2):
        
        num = np.sqrt((GAMMA - 1)*np.pow(M2, 2) + 2)
        den = np.sqrt(GAMMA * (2 * np.pow(M2, 2) - 1) + 1)
        
        self.mach = num / den
    
        state = Shock(self.mach2)
        state.propagate_state2(self)
        self._state2 = state
    

    @property
    def rho2_rho1(self):
        
        num = (GAMMA + 1) * np.pow(super().mach, 2)
        den = (GAMMA - 1)*np.pow(super().mach, 2) + 2
        
        return num / den
    
    @rho2_rho1.setter
    def rho2_rho1(self, r2_r1):
        
        num = 2 * r2_r1
        den = (GAMMA + 1) - r2_r1 * (GAMMA - 1)
        
        M = np.sqrt(num / den)
        self.mach = M
    
    
    @property
    def P2_P1(self):
        
        num = 2*GAMMA*np.pow(super().mach, 2) - (GAMMA - 1)
        den = GAMMA + 1
        
        return num / den
    
    @P2_P1.setter
    def P2_P1(self, P2_P1):
        
        num = P2_P1 * (GAMMA + 1) + (GAMMA - 1)
        den = 2 * GAMMA
        
        M = np.sqrt(num / den)
        self.mach = M
    
    @property
    def T2_T1(self):
        return self.P2_P1 / self.rho2_rho1
    
    @T2_T1.setter
    def T2_T1(self, T2_T1):
        if (T2_T1 <= 1): raise ValueError("T2/T1 must be greater than 1 for Normal Shocks")
        
        g = GAMMA
        x = T2_T1
        
        a = ( g**2 * x**2 + 2*g**2 * x + g**2 + 2*g * x**2 - 12*g * x + 2*g + x**2 + 2*x + 1 ) 
        
        # Full numerator and denominator of the Mach expression 
        num = (g + 1) * np.sqrt(a) + (g**2 + 1) * (x + 1) + 2*g*x - 6*g 
        den = 4 * (g - 1) * g 
        
        M = np.sqrt(num / den) 
        self.mach = M
    
    @property
    def P02_P01(self):
        return self.P2_P1 * np.pow(self.T2_T1, -(GAMMA)/(GAMMA - 1))
    
    @P02_P01.setter
    def P02_P01(self, P02_P01):
        pass

    @property
    def P02_P1(self):
        M = super().mach
        
        t1 = np.pow((GAMMA + 1)*M, 2)/(4*GAMMA*np.pow(M, 2) - 2*(GAMMA - 1))
        t2 = ((1 - GAMMA) + 2*GAMMA*np.pow(M,2)) / (GAMMA + 1)
        
        return np.pow(t1, GAMMA/(GAMMA - 1)) * t2
    
    @P02_P1.setter
    def P02_P1(self, P02_P1):
        
        num = 9*np.pow(GAMMA + 1, 2)
        den = 2*(17*GAMMA + 1)        
        
        K = np.pow(num/den, GAMMA / (GAMMA - 1))
        
        
        a_num = 2*GAMMA*(289*GAMMA + 19)
        a_den = 9*(GAMMA + 1)*(17*GAMMA+1)
        a1 = a_num / a_den * K
        
        b1 = (34*GAMMA)/(3*(GAMMA+1)) * K
        
        c1 = (17*GAMMA + 1)/(GAMMA + 1) * K
        
        a = a1 / 2
        b = b1 - 3*a1
        c = 9/2*a1 - 3*b1 + c1
        
        M = (-b + np.sqrt(np.pow(b, 2) - 4*a*(c - P02_P1)))/2/a
        self.mach = M


    def __setattr__(self, name, value):
        
        match name:
            case 'T2':
                self.T = value / self.T2_T1
                return
            
            case 'P2':
                self.P = value / self.P2_P1
                return
            
            case 'rho2':
                self.rho = value / self.rho2_rho1
                return
        
        
        return super().__setattr__(name, value)


    def __getattr__(self, name):
        
        try:
            match name:
                case 'T2': 
                    return self.T2_T1 * self.T
                
                case 'P2': 
                    return self.P2_P1 * self.P
                
                case 'rho2': 
                    return self.rho2_rho1 * self.rho
                
                case 'P02': 
                    return self.P02_P01 * self.P0
                
                case 'T02': 
                    return self.T0
                
                case 'ds':
                    return -R * np.log(self.P02_P01)
                
                case 'T0_T2':
                    return self.T0 / self.T2
                
                case 'P0_P2':
                    return self.P0 / self.P2
                
                case 'rho0_rho2':
                    return self.rho0 / self.rho2
                
                case "state2":
                    return Isen(self.mach2)
                
        except TypeError:
            raise TypeError(f"Invalid Type: please assign values to normal shock")
                
        
        raise AttributeError(f"Normal has no attirubte \'{name}\'")


# ============================================================ 
# Oblique Shock
# ============================================================

def thetaMachBeta(mach1, beta):
    theta = beta - np.atan((5/(np.pow(mach1*np.sin(beta), 2)) + 1) * np.tan(beta)/6)
    return theta

class Oblique(Shock):
    def __init__(self, M, beta, **kwargs):
        
        if isinstance(M, Normal) or isinstance(M, Oblique):
            M = M.state2
        
        if isinstance(M, Isen):
            state = M
            mach = Q_(state.mach).m
        
            super().__init__(mach, **kwargs)
            self.propagate_state(state)
        
        else:
            super().__init__(M, **kwargs)
        
        self._beta = None       
        self._normal = None
        
        self.beta = beta
        
        
        
        
        
        
    
        
    
    @staticmethod
    def BetaMax(mach):
        epsilon = 1e-10
        beta_max = np.pi / 4
        
        for j in range(20):
            d1 =  diff(lambda B: thetaMachBeta(mach, B), beta_max)
            d2 = diff2(lambda B: thetaMachBeta(mach, B), beta_max)
            
            dB = d1 / d2
            beta_max -= dB
            
            if np.abs(dB) < epsilon: break
        
        beta_max = Q_(beta_max, 'rad').to('deg')
        
        return beta_max
    
    
    @property
    def beta_max(self):
        mach = self.mach.to('').m
        beta_max = Oblique.BetaMax(mach)
        
        return beta_max
    
    @property
    def theta_max(self):
        return thetaMachBeta(self.mach, self.beta_max)
    
    @property
    def beta(self):
        return Q_(self._beta).to("deg")
    
    @property
    def theta(self):
        return thetaMachBeta(self.mach, self.beta).to('deg')
    
    
    def __setattr__(self, name, value):
        
        match name:
            case "beta":
                
                if (Q_(value).to('rad').m > np.pi / 2) or \
                    (Q_(value).to('rad') < self.mu):
                    raise ValueError("Shock is detached")
                
                self._beta = value
                self._normal = Normal(self.mach1n)
                
                return
            
            case "theta":
                
                return
        
        return super().__setattr__(name, value)
    
    
    def __getattr__(self, name):
        
        try:
            match name:
                
                case "mach1n":
                    return self.mach * np.sin(self.beta)
                
                case "mach1t":
                    return self.mach * np.cos(self.beta)
                
                case "mach2n":
                    return np.sqrt((5 + np.pow(self.mach1n, 2)) /\
                        (7*np.pow(self.mach1n, 2) - 1))
                        
                case "mach2t":
                    return self.mach1t * np.sqrt(1 / self.T2_T1)
                
                case "mach2":
                    return np.sqrt(np.pow(self.mach2n, 2) + np.pow(self.mach2t, 2))
                
                case "state2":
                    state = Shock(self.mach2)
                    state.propagate_state2(self)
                    return state
                
                
                case "type":
                    if self.mach2 > 1: return "weak"
                    else: return "strong"
                    
                # ===============
                case 'T2': 
                    return self.T2_T1 * self.T
                
                case 'P2': 
                    return self.P2_P1 * self.P
                
                case 'rho2': 
                    return self.rho2_rho1 * self.rho
                
                case 'P02': 
                    return self.P02_P01 * self.P0
                
                case 'T02': 
                    return self.T0
                
                case 'ds':
                    return -R * np.log(self.P02_P01)
                
                case 'T0_T2':
                    return self.T0 / self.T2
                
                case 'P0_P2':
                    return self.P0 / self.P2
                
                case 'rho0_rho2':
                    return self.rho0 / self.rho2
            
            
            if hasattr(self._normal, name):
                return getattr(self._normal, name)
            
        except TypeError:
            raise TypeError(f"Invalid Type: please assign values to oblique shock")
            
        raise AttributeError(f"Oblique has no attirubte \'{name}\'")
    
    
    @staticmethod
    def IsenBeta(state1: Isen, beta):
        shock1 = Oblique(state1, beta)
        return shock1
    
    @staticmethod
    def IsenTheta(state1: Isen, theta, strong = True):
        
        if isinstance(state1, Shock):
            mach1 = Q_(state1.mach2).m
        else:
            mach1 = Q_(state1.mach).m
        
        
        
        epsilon = 1e-10
        h = 1e-20
        
        mu = np.asin(1 / mach1)
        theta_target = Q_(theta).to('rad').m
        
        if strong: beta1 = mu + h
        else: beta1 = np.pi/2 - h
        
        # beta1 = mu + h
        # beta2 = np.pi / 2 - h
        
        beta_max = Oblique.BetaMax(mach1)
        theta_max = thetaMachBeta(mach1, beta_max)
        
        if theta > theta_max:
            raise ValueError(f"theta exceeds theta max ({theta:.4g} > {theta_max:.4g})")
        
        for j in range(20):
            
            t = thetaMachBeta(mach1, beta1)
            d1 = diff(lambda B: thetaMachBeta(mach1, B), beta1)
            
            dB = (t - theta_target) / d1
            beta1 -= dB
            
            if np.abs(dB) < epsilon: break
            
        # for j in range(20):
            
        #     t = thetaMachBeta(mach1, beta2)
        #     d1 = diff(lambda B: thetaMachBeta(mach1, B), beta2)
            
        #     dB = (t - theta_target) / d1
        #     beta2 -= dB
            
        #     if np.abs(dB) < epsilon: break
            
        theta1 = thetaMachBeta(mach1, beta1)
        
        if np.abs(theta1 - theta_target) > epsilon:
            raise TimeoutError("Could not converge (raise iter range)")
        
        return Oblique(state1, beta1)
        
        
        # theta2 = thetaMachBeta(mach1, beta2)

# ============================================================ 
# Multi-Shock
# ============================================================

class MultiShock:
    
    _variables_A = ("T", "P", "rho")
    _variables_B = ("mu", "beta", "theta", "beta_max", "theta_max")
    
    def __init__(self, state1):
        
        self.state1 = state1
        self.states = [state1]
        
        
    def addShock(self, func, theta = None, beta = None):
        
        state0 = self.states[-1]        
        
        if func is Normal:
            stateN = Normal(state0)
            
        elif func is Oblique:
            
            if theta is not None:
                stateN = Oblique.IsenTheta(state0, theta)
            
            elif beta is not None:
                stateN = Oblique.IsenBeta(state0, beta)
            
            else:
                raise ValueError("Not enough valid inputs (addShock) - 2")
        else:
            print(isinstance(func, Normal.__class__))
            raise ValueError("Not enough valid inputs (addShock) - 1")
        
        self.states.append(stateN)
        
    def __getattr__(self, name):
        
        
        m = re.match(rf"^(mach|{'|'.join([f'{i}0' for i in self._variables_A])}|{'|'.join(self._variables_A)})([0-9]*)$", name)
        
        if m:
            var, i = m.groups()
            i = (int(i)) if (i != '') else 0
            
            if (i == len(self.states)):
                i -= 1
                var = f"{var}2"
            
            return getattr(self.states[i], var)
        
        m = re.match(rf"^({'|'.join(self._variables_B)})([0-9]*)$", name)
        
        if m:
            var, i = m.groups()
            i = (int(i)) if (i != '') else 0
            
            return getattr(self.states[i], var)
        
        m = re.match(rf"^(P0|{'|'.join(self._variables_A)})(\d+)_\1(\d+)$", name)
        
        if m:
            var, i, j = m.groups()
            
            if i == '0':
                return getattr(self.states[int(j)], f"{var}0_{var}")
        
            return getattr(self.states[int(j)], f"{var}2_{var}1")
                
        raise AttributeError(f"unable to access attribute ({name})")