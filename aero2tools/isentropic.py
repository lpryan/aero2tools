from . import core
from .core import *

import numpy as np
import re
 
# ==================== 
# Isentropic Flow 
# ====================

def nuMach(mach):
    sqrtM = np.sqrt(mach**2 - 1)
    nu = np.sqrt(6)*np.atan(sqrtM / np.sqrt(6)) - np.atan(sqrtM)
    return nu



class Isen:
    
    _variables = (
        "mu", "nu", "vel",
        "Tstar_T", "Pstar_P", "rstar_r", "mach_star"
    )
    
    
    def __init__(self, mach: float, **kwargs):
        
        self.GAMMA = kwargs.get("GAMMA", config.GAMMA)
        
        # "hidden variables"
        setattr(self, "_parent", None)
        setattr(self, "_num", None)
        setattr(self, "_updating", False)
        setattr(self, "_generating", True)
            
        # star variables            
        self.Tstar_T = None
        self.Pstar_P = None
        self.rstar_r = None
        self.mach_star = None
        
        self.mach = config.Q_(mach)
        
        # variables
        self.T  = kwargs.get("T",  None)
        self.T0 = kwargs.get("T0", None) if (self.T is None) else None
        
        self.P  = kwargs.get("P",  None)
        self.P0 = kwargs.get("P0", None) if (self.P is None) else None
        
        if (self.T is None) or (self.P is None):
            self.r  = kwargs.get("r",  None)
            self.r0 = kwargs.get("r0", None) if (self.r is None) else None
        else:
            self.r = IdealGas.dens(self.P, self.T)
        
        self.vel = kwargs.get("vel", None)
        
        setattr(self, "_generating", False)
        
        
    
    
    @property
    def mach(self): return self._mach
        
    @mach.setter
    def mach(self, mach: float):
        setattr(self, "_mach", config.Q_(mach))
        
        T0_T = 1 + np.pow(self._mach, 2) * (self.GAMMA - 1)/2
        P0_P = np.pow(T0_T, (self.GAMMA)/(self.GAMMA - 1))
        r0_r = np.pow(T0_T, 1/(self.GAMMA - 1))
        
        k = (self.GAMMA + 1) / 2
        
        Tstar_T = T0_T / k
        Pstar_P = P0_P / np.pow(k, (self.GAMMA)/(self.GAMMA - 1))
        rstar_r = r0_r / np.pow(k, 1/(self.GAMMA - 1))
        mach_star = self._mach / np.sqrt(Tstar_T)
        
        setattr(self, "_T0_T", T0_T)
        setattr(self, "_P0_P", P0_P)
        setattr(self, "_r0_r", r0_r)
        
        setattr(self, "Tstar_T", Tstar_T)
        setattr(self, "Pstar_P", Pstar_P)
        setattr(self, "rstar_r", rstar_r)
        setattr(self, "mach_star", mach_star)
        
        
        num = getattr(self, '_num', None)       
        
        if num is not None:
            parent = getattr(self, "_parent")
            getattr(parent, f"state{3 - num}").reset()
 
    
    # ====================
    # properties    
    #-----
    @property
    def T0_T(self):
        return self._T0_T
    
    @T0_T.setter
    def T0_T(self, T0T: float):
        self.mach = np.sqrt(2/(self.GAMMA - 1) * (T0T - 1))
        
    #-----
    @property
    def P0_P(self):
        return self._P0_P
    
    @P0_P.setter
    def P0_P(self, P0P: float):
        self.mach = np.sqrt(2/(self.GAMMA - 1) * (np.pow(P0P, (self.GAMMA - 1)/(self.GAMMA)) - 1))
    
    #-----
    @property
    def r0_r(self):
        return self._r0_r
    
    @r0_r.setter
    def r0_r(self, r0r: float):
        self.mach = np.sqrt(2/(self.GAMMA - 1) * (np.pow(r0r, (self.GAMMA - 1)) - 1))


    # ====================
    # constructors
    @staticmethod
    def from_temp(T0_T: float | None = None, T0: float | None = None, T: float | None = None, GAMMA: Config | float = config):
        if GAMMA is config: GAMMA = config.GAMMA
        
        if T0_T is None and (T0 is None or T is None): 
            raise ValueError("Not enough inputs (Isen.from_temp)")
        
        elif T0_T is None:
            T0_T = T0 / T
        
        mach = np.sqrt(2 / (GAMMA - 1) * (T0_T - 1))
        
        return Isen(mach, T = T, T0 = T0, GAMMA = GAMMA)
    
    @staticmethod
    def from_pres(P0_P: float | None = None, P0: float | None = None, P: float | None = None, GAMMA: Config | float = config):
        if GAMMA is config: GAMMA = config.GAMMA
        
        if P0_P is None and (P0 is None or P is None):
            raise ValueError("Not enough inputs (Isen.from_pres)")
        
        elif P0_P is None:
            P0_P = P0 / P
            
        mach = np.sqrt(2/(GAMMA - 1) * (np.pow(P0_P, (GAMMA - 1)/(GAMMA)) - 1))
        
        return Isen(mach, P = P, P0 = P, GAMMA = GAMMA)
    
    @staticmethod
    def from_dens(r0_r: float | None = None, r0: float | None = None, r: float | None = None, GAMMA: Config | float = config):
        if GAMMA is config: GAMMA = config.GAMMA
        
        if r0_r is None and (r0 is None or r is None):
            raise ValueError("Not enough inputs (Isen.from_dens)")
        
        elif r0_r is None:
            r0_r = r0 / r
            
        mach = np.sqrt(2/(GAMMA - 1) * (np.pow(r0_r, (GAMMA - 1)/(GAMMA)) - 1))
        
        return Isen(mach, r = r, r0 = r0, GAMMA = GAMMA)
    
    @staticmethod
    def from_vel(temp: float, vel: float, GAMMA: Config | float = config):
        if GAMMA is config: GAMMA = config.GAMMA
        
        mach = SpeedOfSound.mach(temp, vel)
        
        return Isen(mach, T = temp, GAMMA = GAMMA)

    
    # ====================
    # attribute modifier/accessor
    
    def __setattr__(self, name, value):
        updating = False
        
        
        if name.startswith("_") or name in ("GAMMA",):
            super().__setattr__(name, value)
            return
        
        if name in self._variables and getattr(self, "_generating", False):
            super().__setattr__(name, value)
            return
        
        
        if name in ("T", "T0", "P", "P0", "r", "r0"):
            
            if getattr(self, "_generating", False) and (getattr(self, name, None) is not None): return            
            
            super().__setattr__(name, value)
            
            if value is None: return
            
            setattr(self, "_updating", True)
            
            try:
                
                num = getattr(self, '_num', None)
        
                if num is not None:
                    parent = getattr(self, "_parent")
                    other = getattr(parent, f"state{3 - num}")
                    
                    updating = not getattr(other, "_updating", False)
                
                match name:
                    # Temperature
                    case "T":
                        super().__setattr__("T0", value * self.T0_T)
                        super().__setattr__("Tstar", value * self.Tstar_T)
                        
                        if updating: other.T = value * np.pow(getattr(parent, "T2_T1"), 3 - 2*num)
                        
                    case "T0":
                        self.T = value / self.T0_T
                    
                    # Pressure
                    case "P":
                        super().__setattr__("P0", value * self.P0_P)
                        super().__setattr__("Pstar", value * self.Pstar_P)
                        
                        if updating: other.P = value * np.pow(getattr(parent, "P2_P1"), 3 - 2*num)
                        
                    case "P0":
                        self.P = value / self.P0_P
                        
                    # Density
                    case "r":
                        super().__setattr__("r0", value * self.r0_r)
                        super().__setattr__("rstar", value * self.rstar_r)
                        
                        if updating: other.r = value * np.pow(getattr(parent, "r2_r1"), 3 - 2*num)
                    
                    case "r0":
                        self.r = value / self.r0_r
                    
            finally:
                setattr(self, "_updating", False)
        
        
        if getattr(self, "T", None) is not None:
            super().__setattr__("vel", SpeedOfSound.vel(self.T, self.mach))
                    
        validIdealGas = [i for i in ("T", "P", "r") if getattr(self, f"{i}", None) is None]
        
        if validIdealGas.count(None) == 1:    
            match validIdealGas[0]:
                case "T":
                    self.T = IdealGas.temp(self.P, self.r)
                case "P":
                    self.P = IdealGas.pres(self.T, self.r)
                case "r":
                    self.r = IdealGas.dens(self.P, self.T)
            return
        
        if name in dir(self):
            super().__setattr__(name, value)
            return
        
        raise AttributeError(f"Invalid assignment ({name})")
        
    def __getattr__(self, name):
        
        if name not in ["mach", "_mach"]: M = config.Q_(self.mach).to('').m
        
        match name:
            
            case "nu": 
                if M >= 1:
                    nu = nuMach(self.mach)
                    return nu.to('deg')
            
            case "mu": 
                if M >= 1:
                    return np.asin(1 / self.mach)
        
        raise AttributeError("Isen does not have this attribute")
    
    
    def reset(self):
        
        num = getattr(self, '_num', None)
        parent = getattr(self, "_parent")
        
        other = getattr(parent, f"state{3 - num}")
        
        if other.T is not None:
            self.T = other.T * np.pow(getattr(parent, "T2_T1"), 2*num - 3)
        
        if other.P is not None:
            self.P = other.P * np.pow(getattr(parent, "P2_P1"), 2*num - 3)
        
        if self.T and self.P is not None:
            self.r = IdealGas.dens(self.P, self.T)
        
        
    
    
    def propagate_state(self, state):
        
        self.mach = state.mach
        
        try: self.T = state.T
        except TypeError: pass
        
        try: self.P = state.P
        except TypeError: pass
        
        try: self.r = state.r
        except TypeError: pass

    
    # ====================
    # string
    def __str__(self):
        return "-"*3+"Isentropic"+"-"*3+f"\nmach: {self.mach:>6.3f~P}\tT0/T: {self.T0_T:>6.3f~P}" + \
               f"\nP0/P: {self.P0_P:>6.3f~P}\tρ0/ρ: {self.r0_r:>6.3f~P}"



# ============================================================ 
# Isentropic Transfer Class
# - transfer between any two isentropic states
# ============================================================

class IsenTranslate:
    
    
    def __init__(self, state1: Isen, state2: Isen):
        
        self.state1 = state1
        self.state2 = state2
        
        setattr(self.state1, "_parent", self)
        setattr(self.state1, "_num", 1)
        
        setattr(self.state2, "_parent", self)
        setattr(self.state2, "_num", 2)
        
        self.state1.reset()
        self.state2.reset()
    
    # ================
    #   Constructors
    @staticmethod
    def from_temp1(state1: Isen, T2_T1 = None, T1 = None, T2 = None):
        
        if T1 is not None: state1.T = T1
        T1 = getattr(state1, "T", T1)
        
        if (T2_T1 is None) and (T1 is None or T2 is None):
            raise ValueError("Not enough inputs (IsenTranslate.from_temp1)")
        
        elif not (T1 is None or T2 is None):
            T2_T1 = T2 / T1
        
        T0_T2 = state1.T0_T / T2_T1
        state2 = Isen.from_temp(T0_T2, T = T2)
        
        inter = IsenTranslate(state1, state2)
        inter.state2.reset()
            
        return inter
    
    @staticmethod
    def from_temp2(state2: Isen, T2_T1 = None, T1 = None, T2 = None):
        
        if T2 is not None: state2.T = T2
        T2 = getattr(state1, "T", T2)
        
        if (T2_T1 is None) and (T1 is None or T2 is None):
            raise ValueError("Not enough inputs (IsenTranslate.from_temp2)")
        
        elif not (T1 is None or T2 is None):
            T2_T1 = T2 / T1
            
        T0_T1 = state2.T0_T * T2_T1
        state1 = Isen.from_temp(T0_T1, T = T1)
        
        inter = IsenTranslate(state1, state2)
        inter.state1.reset()
        
        return inter
    
    @staticmethod
    def from_pres1(state1: Isen, P2_P1 = None, P1 = None, P2 = None):
        
        if P1 is not None: state1.P = P1
        P1 = getattr(state1, "P", P1)
        
        if (P2_P1 is None) and (P1 is None or P2 is None):
            raise ValueError("Not enough inputs (IsenTranslate.from_pres1)")
        
        elif not (P1 is None or P2 is None):
            P2_P1 = P2 / P1
            
        P0_P2 = state1.P0_P / P2_P1
        state2 = Isen.from_pres(P0_P2, P = P2)
        
        inter = IsenTranslate(state1, state2)
        inter.state2.reset()
        
        return inter
    
    @staticmethod
    def from_pres2(state2: Isen, P2_P1 = None, P1 = None, P2 = None):
        
        if P2 is not None: state2.P = P2
        P2 = getattr(state2, "P", P2)
        
        if (P2_P1 is None) and (P1 is None or P2 is None):
            raise ValueError("Not enough inputs (IsenTranslate.from_pres2)")
        
        elif not (P1 is None or P2 is None):
            P2_P1 = P2 / P1
            
        P0_P1 = state2.P0_P / P2_P1
        state1 = Isen.from_pres(P0_P1, P = P1)
        
        inter = IsenTranslate(state1, state2)
        inter.state1.reset()
        
        return inter
    
    @staticmethod
    def from_dens1(state1: Isen, r2_r1 = None, r1 = None, r2 = None):
        
        if r1 is not None: state1.r = r1
        r1 = getattr(state1, "r", r1)
        
        if (r2_r1 is None) and (r1 is None or r2 is None):
            raise ValueError("Not enough inputs (IsenTranslate.from_dens1)")
        
        elif not (r1 is None or r2 is None):
            r2_r1 = r2 / r1
            
        r0_r2 = state1.r0_r / r2_r1
        state2 = Isen.from_dens(r0_r2, r = r2)
        
        inter = IsenTranslate(state1, state2)
        inter.state2.reset()
        
        return inter
    
    @staticmethod
    def from_dens2(state2: Isen, r2_r1 = None, r1 = None, r2 = None):
        
        if r2 is not None: state2.r = r2
        r2 = getattr(state2, "r", r2)
        
        if (r2_r1 is None) and (r1 is None or r2 is None):
            raise ValueError("Not enough inputs (IsenTranslate.from_dens2)")
        
        elif not (r1 is None or r2 is None):
            r2_r1 = r2 / r1
            
        r0_r1 = state2.r0_r / r2_r1
        state1 = Isen.from_pres(r0_r1, r = r1)
        
        inter = IsenTranslate(state1, state2)
        inter.state1.reset()
        
        return inter
    
    @staticmethod
    def from_dtheta1(state1: Isen | float, dtheta, **kwargs):
        
        if isinstance(state1, float):
            state1 = Isen(state1, **kwargs)

        nu1 = state1.nu
        nu2 = dtheta + nu1
        
        mach2 = optimize.target(nuMach, 2, nu2)
        state2 = Isen(mach2)
        
        exp = IsenTranslate(state1, state2)
        
        return exp
    
    
    # ====================
    # attribute accessor
    
    @property
    def T2_T1(self):
        return self.state1.T0_T / self.state2.T0_T
        
    @property
    def P2_P1(self):
        return self.state1.P0_P / self.state2.P0_P
    
    @property
    def r2_r1(self):
        return self.state1.r0_r / self.state2.r0_r
    
    @property
    def P02_P01(self):
        return self.P2_P1 * np.pow(1 / self.T2_T1, 7/2)
    
    
    def __getattr__(self, name):
        
        m = re.match(rf"^(P|T|r)0_\1(1|2)$", name)
        
        if m:
            var, i = m.groups()
            
            match var:
                case "P": 
                    return getattr(self, f"state{i}").P0_P
                case "T": 
                    return getattr(self, f"state{i}").T0_T
                case "r": 
                    return getattr(self, f"state{i}").r0_r
                
        
        m = re.match(rf"^((P|T|r)(star|0)?|mach|vel)(|1|2)$", name)
        
        if m:
            var, i, j, k = m.groups()
            if k == '': k = 1
            return getattr(getattr(self, f"state{k}"), f"{var}", None)
        
        if name == "entropy":
            return config.CP * np.log(self.T2_T1) - config.R * np.log(self.P2_P1)
        
        if name == "dtheta":
            return self.state2.nu - self.state1.nu
        
        raise AttributeError(f"{name} not found")