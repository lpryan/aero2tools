from . import core
from .core import *
import numpy as np

# ============================================================ 
# Isentropic Flow 
# ============================================================

class Isen:
    
    _variables = ("vel", "mu", "nu", 
                  "Tstar_T", "Pstar_P", "rhostar_rho", "mach_star",
                  "Tstar", "Pstar", "rhostar")
    
    def __init__(self, M: float, **kwargs):
        
        self.GAMMA = kwargs.get("GAMMA")
        if kwargs.get("GAMMA") is None: self.GAMMA = core.GAMMA
        
        self._parent = None
        self._num = None
        self._updating = False
        
        self._T0_T = None
        self._P0_P = None
        self._rho0_rho = None
        
        self.T = None
        self.T0 = None
        self.P = None
        self.P0 = None
        self.rho = None
        self.rho0 = None
        
        self.vel = None
        self.mu = None
        self.nu = None
        
        self.Tstar_T = None
        self.Pstar_P = None
        self.rhostar_rho = None
        self.mach_star = None
        
        self.Tstar = None
        self.Pstar = None
        self.rhostar = None
        
        # Initialize
        self.mach = Q_(M)
    

    # =========================
    # Property Management
    @property
    def mach(self):
        return self._mach
    
    @mach.setter
    def mach(self, M):
        
        setattr(self, "_mach", Q_(M))
        
        T0_T = 1 + np.pow(self._mach, 2) * (self.GAMMA - 1)/2
        P0_P = np.pow(T0_T, (self.GAMMA)/(self.GAMMA - 1))
        rho0_rho = np.pow(T0_T, 1/(self.GAMMA - 1))
        
        k = (self.GAMMA + 1)/2
        
        Tstar_T = T0_T / k
        Pstar_P = P0_P / np.pow(k, (self.GAMMA)/(self.GAMMA - 1))
        rhostar_rho = rho0_rho / np.pow(k, 1/(self.GAMMA - 1))
        mach_star = self._mach / np.sqrt(Tstar_T)
        
        setattr(self, "_T0_T", T0_T)
        setattr(self, "_P0_P", P0_P)
        setattr(self, "_rho0_rho", rho0_rho)
        
        setattr(self, "Tstar_T", Tstar_T)
        setattr(self, "Pstar_P", Pstar_P)
        setattr(self, "rhostar_rho", rhostar_rho)
        setattr(self, "mach_star", mach_star)
        
        
        if self._mach >= 1:
            self.mu = np.asin(1 / self._mach)
            self.nu = Q_(20, 'deg')        
        
        if self._num is not None:
            getattr(self._parent, f"state{3 - self._num}").reset()
        
        
    
    #-----
    @property
    def T0_T(self):
        return self._T0_T
    
    @T0_T.setter
    def T0_T(self, T0T):
        self.mach = np.sqrt(2/(self.GAMMA - 1) * (T0T - 1))
        
    #-----
    @property
    def P0_P(self):
        return self._P0_P
    
    @P0_P.setter
    def P0_P(self, P0P):
        self.mach = np.sqrt(2/(self.GAMMA - 1) * (np.pow(P0P, (self.GAMMA - 1)/(self.GAMMA)) - 1))
    
    #-----
    @property
    def rho0_rho(self):
        return self._rho0_rho
    
    @rho0_rho.setter
    def rho0_rho(self, rho0rho):
        self.mach = np.sqrt(2/(GAMMA - 1) * (np.pow(rho0rho, (GAMMA - 1)) - 1))
        
    # =========================
    # Constructors
    @staticmethod
    def from_temp(T0_T = None, T0 = None, T = None, GAMMA = core.GAMMA):
        
        if T0_T is None and (T0 is None or T is None):
            print("Invalid Input (No information for Isen.from_temp)")
            return -1
        
        elif T0_T is None:
            T0_T = T0 / T

        M = np.sqrt(2/(GAMMA - 1) * (T0_T - 1))
        return Isen(M)
    
    @staticmethod
    def from_pressure(P0_P = None, P0 = None, P = None, GAMMA = core.GAMMA):

        if P0_P is None and (P0 is None or P is None):
            print("Invalid Input (No information for Isen.from_pressure)")
            return -1
        
        elif P0_P is None:
            P0_P = P0 / P

        M = np.sqrt(2/(GAMMA - 1) * (np.pow(P0_P, (GAMMA - 1)/(GAMMA)) - 1))
        return Isen(M)
    
    @staticmethod
    def from_density(rho0_rho = None, rho0 = None, rho = None, GAMMA = core.GAMMA):

        if rho0_rho is None and (rho0 is None or rho is None):
            print("Invalid Input (No information for Isen.from_pressure)")
            return -1
        
        elif rho0_rho is None:
            rho0_rho = rho0 / rho

        M = np.sqrt(2/(GAMMA - 1) * (np.pow(rho0_rho, (GAMMA - 1)) - 1))
        return Isen(M)
    
    @staticmethod
    def from_velocity(temp, velocity, GAMMA = core.GAMMA):
        M = SpeedOfSound.mach(temp, velocity, GAMMA)
        
        state = Isen(M)
        state.T = temp
        
        return state        
    
    
    # =========================
    # Attribute Setter
    def __setattr__(self, name, value):
        updating = False
        
        # --- allow internal / private attributes to be set normally ---
        if name.startswith("_") or name in ("GAMMA",):
            super().__setattr__(name, value)
            return
        
        if name in self._variables:
            super().__setattr__(name, value)
            return
                    
        # --- intercept state assignments ---
        if name in ("T", "T0", "P", "P0", "rho", "rho0"):
            super().__setattr__(name, value)
            
            if value is None: return
            
            self._updating = True
            try:
                
                if self._num is not None:
                    other = getattr(self._parent, f"state{3 - self._num}")
                    updating = (self._num is not None) and (not getattr(other, "_updating", False))
                
                match name:
                    # Temperature
                    case "T":
                        super().__setattr__("T0", value * self._T0_T)
                        super().__setattr__("Tstar", value * self.Tstar_T)
                        
                        if updating: other.T = value * np.pow(getattr(self._parent, "T2_T1"), 3 - 2*self._num)
                                            
                    case "T0":
                        self.__setattr__("T", value / self._T0_T)
                    
                    # Pressure
                    case "P":
                        super().__setattr__("P0", value * self._P0_P)
                        super().__setattr__("Pstar", value * self.Pstar_P)
                        
                        if updating: other.P = value * np.pow(getattr(self._parent, "P2_P1"), 3 - 2*self._num)
                        
                    case "P0":
                        self.__setattr__("P", value / self._P0_P)
                    
                    # Density
                    case "rho":
                        super().__setattr__("rho0", value * self._rho0_rho)
                        super().__setattr__("rhostar", value * self.rhostar_rho)
                        
                        if updating: other.rho = value * np.pow(getattr(self._parent, "rho2_rho1"), 3 - 2*self._num)
                        
                    case "rho0":
                        self.__setattr__("rho", value / self._rho0_rho)
        
            finally:
                self._updating = False
        
        
        if None not in [self.T, self.P] and self.rho is None:
            self.rho = IdealGas.density(self.P, self.T)
            return
        
        if None not in [self.T, self.rho] and self.P is None:
            self.P = IdealGas.pressure(self.T, self.rho)
            return
            
        if None not in [self.rho, self.P] and self.T is None:
            self.T = IdealGas.temperature(self.P, self.rho)
            return
        
        if self.T is not None:
            super().__setattr__("vel", SpeedOfSound.velocity(self.T, self.mach))
        
        
        if name in dir(self):
            super().__setattr__(name, value)
            return
        
        
        # Default
        raise AttributeError(f"Invalid assignment ({name})")
    
    
    def reset(self):
        
        other = getattr(self._parent, f"state{3 - self._num}")
        
        self.T = other.T * np.pow(getattr(self._parent, "T2_T1"), 2*self._num - 3)
        self.P = other.P * np.pow(getattr(self._parent, "P2_P1"), 2*self._num - 3)
        
        pass
    
    def propagate_state(self, state):
    
        try: self.T = state.T
        except TypeError: pass
        
        try: self.P = state.P
        except TypeError: pass
        
        try: self.rho = state.rho
        except TypeError: pass
    
    
    # =========================
    # String
    def __str__(self):
        return "-"*3+"Isentropic"+"-"*3+f"\nmach: {self.mach:>6.3f~P}\tT0/T: {self.T0_T:>6.3f~P}" + \
               f"\nP0/P: {self.P0_P:>6.3f~P}\tρ0/ρ: {self.rho0_rho:>6.3f~P}"






# ============================================================ 
# Isentropic Transfer Class
# - transfer between any two isentropic states
# ============================================================

class IsenTranslate:
    
    ratio_names = ("T2_T1", "P2_P1", "rho2_rho1", "P02_P01")
    
    def __init__(self, state1: Isen, state2: Isen):
        
        self.state1 = state1
        self.state2 = state2
        
        state1._parent = self
        state1._num = 1
        
        state2._parent = self
        state2._num = 2
        
        self._T2_T1 = None
        self._P2_P1 = None
        self._rho2_rho1 = None
        self._P02_P01 = None
        
    # ======================
    #      Constructors
    @staticmethod
    def from_temp1(state1: Isen, T2_T1 = None, T1 = None, T2 = None):
        
        if state1.T is not None:
            T1 = state1.T
        
        if (T2_T1 is None) and (T1 is None or T2 is None):
            raise ValueError("Not enough inputs (from_temp1)")
        
        elif not (T1 is None or T2 is None):
            T2_T1 = T2 / T1
        
        T0_T2 = state1.T0_T / T2_T1
        state2 = Isen.from_temp(T0_T2)
        
        inter = IsenTranslate(state1, state2)
        
        
        if not (T1 is None or T2 is None):
            inter.state1.T = T1
            
        return inter
    
    @staticmethod
    def from_temp2(state2: Isen, T2_T1 = None, T1 = None, T2 = None):
        
        if state2.T is not None:
            T2 = state2.T
            
        if (T2_T1 is None) and (T1 is None or T2 is None):
            raise ValueError("Not enough inputs (from_temp1)")
        
        elif not (T1 is None or T2 is None):
            T2_T1 = T2 / T1
            
        T0_T1 = state2.T0_T * T2_T1
        state1 = Isen.from_temp(T0_T1)
        
        inter = IsenTranslate(state1, state2)
        
        if not (T1 is None or T2 is None):
            inter.state1.T = T1
            
        return inter
    
    @staticmethod
    def from_pressure1(state1: Isen, P2_P1 = None, P1 = None, P2 = None):
        
        if state1.P is not None:
            P1 = state1.P
        
        if (P2_P1 is None) and (P1 is None or P2 is None):
            raise ValueError("Not enough inputs (from_pressure1)")
        
        elif not (P1 is None or P2 is None):
            P2_P1 = P2 / P1
        
        P0_P2 = state1.P0_P / P2_P1
        state2 = Isen.from_pressure(P0_P2)
        
        inter = IsenTranslate(state1, state2)
        
        
        if not (P1 is None or P2 is None):
            inter.state1.P = P1
            
        return inter
    
    @staticmethod
    def from_pressure2(state2: Isen, P2_P1 = None, P1 = None, P2 = None):
        
        if state2.P is not None:
            P2 = state2.P
            
        if (P2_P1 is None) and (P1 is None or P2 is None):
            raise ValueError("Not enough inputs (from_pressure2)")
        
        elif not (P1 is None or P2 is None):
            P2_P1 = P2 / P1
            
        P0_P1 = state2.P0_P * P2_P1
        state1 = Isen.from_pressure(P0_P1)
        
        inter = IsenTranslate(state1, state2)
        
        if not (P1 is None or P2 is None):
            inter.state1.P = P1
            
        return inter
    
    
    @staticmethod 
    def from_density1(state1: Isen, rho2_rho1 = None, rho1 = None, rho2 = None): 
        
        if state1.rho is not None: 
            rho1 = state1.rho 
            
        if (rho2_rho1 is None) and (rho1 is None or rho2 is None): 
            raise ValueError("Not enough inputs (from_temp1)") 
        
        elif not (rho1 is None or rho2 is None): 
            rho2_rho1 = rho2 / rho1 
            
        rho0_rho2 = state1.rho0_rho / rho2_rho1 
        state2 = Isen.from_temp(rho0_rho2) 
        
        inter = IsenTranslate(state1, state2) 
        
        if not (rho1 is None or rho2 is None): 
            inter.state1.rho = rho1
        
        return inter
    
    
    @staticmethod
    def from_density2(state2: Isen, rho2_rho1 = None, rho1 = None, rho2 = None):
        
        if state2.rho is not None:
            rho2 = state2.rho
            
        if (rho2_rho1 is None) and (rho1 is None or rho2 is None):
            raise ValueError("Not enough inputs (from_temp1)")
        
        elif not (rho1 is None or rho2 is None):
            rho2_rho1 = rho2 / rho1
            
        rho0_rho1 = state2.rho0_rho * rho2_rho1
        state1 = Isen.from_temp(rho0_rho1)
        
        inter = IsenTranslate(state1, state2)
        
        if not (rho1 is None or rho2 is None):
            inter.state1.rho = rho1
        
        return inter
    
            
    def __getattr__(self, name):
        
        if name in IsenTranslate.ratio_names:
            self.update()
            return getattr(self, f"_{name}")
        
        if name == "entropy":
            self.update()
            return CP * np.log(self._T2_T1) - R * np.log(self._P2_P1)
        
        raise AttributeError(f"{name} not found")
    
    def update(self):
        
        T2_T1 = self.state1.T0_T / self.state2.T0_T
        P2_P1 = self.state1.P0_P / self.state2.P0_P
        rho2_rho1 = self.state1.rho0_rho / self.state2.rho0_rho
        P02_P01 = P2_P1 * np.pow(1 / T2_T1, 7/2)
        
        self._T2_T1 = T2_T1
        self._P2_P1 = P2_P1
        self._rho2_rho1 = rho2_rho1
        self._P02_P01 = P02_P01