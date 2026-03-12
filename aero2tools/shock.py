from .core import *
from .isentropic import *

import numpy as np
import re


# ============================================================ 
# General Shock
# ============================================================

class Shock(Isen):
    
    def __init__(self, M, **kwargs):
        
        if isinstance(M, Normal) or isinstance(M, Oblique):
            M = M.state2
            
        if isinstance(M, Isen):
            state = M
            mach = config.Q_(state.mach).m
            
            if (mach < 1): raise ValueError(f"shock cannot occure with subsonic mach (M = {M})")
            
            super().__init__(mach, **kwargs)
            self.propagate_state(state)
            
        else:
            if (M < 1): raise ValueError(f"shock cannot occure with subsonic mach (M = {M})")
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
            raise TypeError(f"Invalid Type: please assign values to shock")
        
        try:
            return super().__getattr__(name)
        
        except AttributeError:
            raise AttributeError(f"Shock has no attirubte \'{name}\'")

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
        
    def __getattr__(self, name):
        
        try:
            return super().__getattr__(name)
        
        except AttributeError:
            raise AttributeError(f"Normal has no attirubte \'{name}\'")
        

# ============================================================ 
# Oblique Shock
# ============================================================

def thetaMachBeta(mach1, beta):
    theta = beta - np.atan((5/(np.pow(mach1*np.sin(beta), 2)) + 1) * np.tan(beta)/6)
    return theta

class Oblique(Shock):
    
    def __init__(self, M, beta, **kwargs):
        super().__init__(M, **kwargs)
        
        self._normal = None
        self.beta = beta
    
    
    def BetaMax(mach = None):
        mach = config.Q_(mach).to('').m
        beta_max = optimize.optimize(lambda B: thetaMachBeta(mach, B), np.pi / 4)
        return config.Q_(beta_max,'rad')
    
    @property
    def beta_max(self):
        return Oblique.BetaMax(self.mach)
    
    @property
    def theta_max(self):
        return thetaMachBeta(self.mach, self.beta_max)
    
    @property
    def beta(self):
        return config.Q_(getattr(self, "_beta", None)).to('deg')
    
    @beta.setter
    def beta(self, B):
        
        angle = config.Q_(B).to('rad').m
        
        if (angle > np.pi / 2) or (angle < self.mu.to('rad').m):
            raise ValueError("Shock is detached")
        
        setattr(self, "_beta", B)
        setattr(self, "_normal", Normal(self.mach1n))
    
    @property
    def theta(self):
        return thetaMachBeta(self.mach, self.beta).to('deg')
    
    
    @property
    def T2_T1(self):
        return self._normal.T2_T1

    @property
    def P2_P1(self):
        return self._normal.P2_P1

    @property
    def r2_r1(self):
        return self._normal.r2_r1

    @property
    def P02_P01(self):
        return self._normal.P02_P01
    
    
    def __getattr__(self, name):
        
        try:
            match name:
                
                case "mach1n":
                    return self.mach * np.sin(self.beta)
                
                case "mach1t":
                    return self.mach * np.cos(self.beta)
                
                case "mach2n":
                    return self._normal.mach2
                        
                case "mach2t":
                    return self.mach1t * np.sqrt(1 / self.T2_T1)
                
                case "mach2":
                    return np.sqrt(self.mach2n**2 + self.mach2t**2)
                
                case "state2":
                    state = Shock(self.mach2)
                    state.propagate_state2(self)
                    return state
                
                case "type":
                    return "weak" if self.mach2 < 1 else "strong"
        
            return super().__getattr__(name)
        
        except TypeError:
            raise TypeError(f"Invalid Type: please assign values to oblique shock")
        
        except AttributeError:
            raise AttributeError(f"Oblique has no attirubte \'{name}\'")
        
    
    @staticmethod
    def IsenBeta(state1: Isen, beta, **kwargs):
        shock1 = Oblique(state1, beta, **kwargs)
        return shock1
    
    @staticmethod
    def IsenTheta(state1: Isen | float, theta, strong = True, **kwargs):
        Q_ = config.Q_
        
        if isinstance(state1, Shock):
            mach1 = Q_(state1.mach2).m

        elif isinstance(state1, Isen):
            mach1 = Q_(state1.mach).m
        
        else:
            
            try:
                mach1 = Q_(float(state1)).to('').m
            
            except ValueError:
                raise ValueError(f"invalid input for mach {state1}")
            
            state1 = Isen(mach1, **kwargs)
            
            
        epsilon = 1e-10
        h = 1e-20
        
        mu = np.asin(1 / mach1)
        theta_target = Q_(theta).to('rad').m
        
        beta1 = (mu + h) if strong else (np.pi/2 - h)
        
        beta_max = Oblique.BetaMax(mach = mach1)
        theta_max = thetaMachBeta(mach1, beta_max)
        
        if theta > theta_max:
            raise ValueError(f"theta exceeds maximum ({theta:~.4g} > {theta_max:~.4g})")
        
        beta_opt = optimize.target(lambda B: thetaMachBeta(mach1, B), beta1, theta_target)
        theta_opt = thetaMachBeta(mach1, beta_opt)
        
        if np.abs(theta_opt - theta_target) > epsilon:
            raise TimeoutError("Could not converge (raise iter range)")
        
        shock1 = Oblique(state1, beta_opt)
        return shock1

          



    
# ============================================================ 
# Isen Tracker
# ============================================================  

Expansion = IsenTranslate

class IsenTracker:
    
    _variables_A = ("T", "P", "r") # variables with ratios
    _variables_B = ("mu", "nu", "beta", "theta", "beta_max", "theta_max", "fwd", "rwd") # variables without ratios
    
    def __init__(self, state0):
        
        self.state0 = state0
        self.states = [state0]
    
    def addShock(self, func, theta = None, beta = None, **kwargs):
        
        statePrev = self.states[-1]
        
        if func is Normal:
            stateNew = Normal(statePrev, **kwargs)
            
        elif func is Oblique:
            
            if (beta is not None) and (theta is not None):
                raise ValueError("Too many inputs (addShock: Oblique)")                
                
            elif beta is not None:
                stateNew = Oblique.IsenBeta(statePrev, beta, **kwargs)
            
            elif theta is not None:
                stateNew = Oblique.IsenTheta(statePrev, theta, True, **kwargs)
            
            else:
                raise ValueError("Not enough valid inputs (addShock: Oblique)")
        
        elif func is Expansion:
            
            T2 = kwargs.get('T2')
            P2 = kwargs.get('P2')
            r2 = kwargs.get('r2')
        
            if theta is not None:
                stateNew = Expansion.from_theta1(statePrev, theta)
            
            elif T2 is not None:
                stateNew = Expansion.from_temp1(statePrev, T2 = T2)
            
            elif P2 is not None:
                stateNew = Expansion.from_pres1(statePrev, P2 = P2)
            
            elif r2 is not None:
                stateNew = Expansion.from_dens1(statePrev, r2 = r2)
                
            else:
                raise ValueError("Not enough valid inputs (addShock: Expansion)")
        
        else:
            raise ValueError("Not enough valid inputs (addShock)")
        
        self.states.append(stateNew)
            
    
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
            
            if (i == len(self.states)):
                i -= 1
                var = f"{var}2"
            
            return getattr(self.states[i], var)
        
        m = re.match(rf"^(P0|{'|'.join(self._variables_A)})(\d+)_\1(\d+)$", name)
        
        if m:
            var, i, j = m.groups()
            
            if i == '0':
                return getattr(self.states[int(j)], f"{var}0_{var}")
            
            return getattr(self.states[int(j)], f"{var}2_{var}1")
                
        raise AttributeError(f"unable to access attribute ({name})")
        
        
    
    
    
    
    













