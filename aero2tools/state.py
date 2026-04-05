from __future__ import annotations
from .core import *
from .optimize import *

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .Relations.relation import Relation
    from .Relations.interIsen import InterIsen
    from .Relations.normal import Normal

import numpy as np
import re

global num
num = 0


# ===========================
# Isentropic State
# ===========================

def nuMach(mach, GAMMA = config):
    if isinstance(GAMMA, Config):
        GAMMA = GAMMA.GAMMA
    sqrtM = np.sqrt(mach**2 - 1)
    sqrtGamma = np.sqrt((GAMMA + 1) / (GAMMA - 1))
    nu = sqrtGamma * np.atan(sqrtM / sqrtGamma) - np.atan(sqrtM)
    return nu

def AA0Mach(mach, GAMMA = config):
    if isinstance(GAMMA, Config):
        GAMMA = GAMMA.GAMMA    
    
    k = (GAMMA + 1)/2
    T0T = 1 + (GAMMA - 1)/2 * mach**2
    AA0t = np.pow(T0T/k, k/(GAMMA - 1)) / mach
    return AA0t


class Isen:
    
    def __init__(self, M, **kwargs):
        self.GAMMA = kwargs.get("GAMMA", config.GAMMA)
        setattr(self, 'UUID', uuid.uuid4())

        setattr(self, "_parent", [])
        setattr(self, "_children", [])
        setattr(self, "_tracker", None)
        setattr(self, "_updating", False)
        setattr(self, "_generating", True)
        
        # star variables
        self.Tstar_T = None
        self.Pstar_P = None
        self.rstar_r = None
        self.mach_star = None
        
        # set mach
        self.mach = config.Q_(M)
        
        # variables
        self.T  = kwargs.get("T", None)        
        self.P  = kwargs.get("P",  None)
        
        if (self.T is None) or (self.P is None):
            self.r  = kwargs.get("r",  None)
        else:
            self.r = IdealGas.dens(self.P, self.T)
            
        self.A  = kwargs.get("A",  None)
        self.A0 = kwargs.get("A0", None) if (self.A is None) else None
        
        if getattr(self, 'vel', None) is None:
            self.vel = kwargs.get("vel", None)
        
        setattr(self, "_generating", False)
    
    
    # ----------------------------------
    # Properties
    # ----------------------------------
    
    # -- mach --    
    @property
    def mach(self):
        return self._mach
    
    @mach.setter
    @config.wrap(None, (None, ''), False)
    def mach(self, mach: float):
        setattr(self, "_mach", config.Q_(mach))
        
        T0_T = 1 + np.pow(self._mach, 2) * (self.GAMMA - 1)/2
        P0_P = np.pow(T0_T, (self.GAMMA)/(self.GAMMA - 1))
        r0_r = np.pow(T0_T, 1/(self.GAMMA - 1))
        
        k = (self.GAMMA + 1) / 2
        
        A_A0 = AA0Mach(self._mach, self.GAMMA)
        
        Tstar_T = T0_T / k
        Pstar_P = P0_P / np.pow(k, (self.GAMMA)/(self.GAMMA - 1))
        rstar_r = r0_r / np.pow(k, 1/(self.GAMMA - 1))
        Astar_A = np.pow(Tstar_T/k, k/(self.GAMMA - 1))
        mach_star = self._mach / np.sqrt(Tstar_T)
        
        setattr(self, "_T0_T", T0_T)
        setattr(self, "_P0_P", P0_P)
        setattr(self, "_r0_r", r0_r)
        
        setattr(self, "Tstar_T", Tstar_T)
        setattr(self, "Pstar_P", Pstar_P)
        setattr(self, "rstar_r", rstar_r)
        setattr(self, "mach_star", mach_star)
        
        setattr(self, "_A_A0", A_A0)
        
        self.attr_pulse()
    
    # -- T0/T --
    @property
    def T0_T(self):
        return self._T0_T
    
    @T0_T.setter
    def T0_T(self, T0T: float):
        self.mach = np.sqrt(2/(self.GAMMA - 1) * (T0T - 1))

    # -- P0/P --
    @property
    def P0_P(self):
        return self._P0_P
    
    @P0_P.setter
    def P0_P(self, P0P: float):
        self.mach = np.sqrt(2/(self.GAMMA - 1) * (np.pow(P0P, (self.GAMMA - 1)/(self.GAMMA)) - 1))
        
    # -- r0/r --
    @property
    def r0_r(self):
        return self._r0_r
    
    @r0_r.setter
    def r0_r(self, r0r: float):
        self.mach = np.sqrt(2/(self.GAMMA - 1) * (np.pow(r0r, (self.GAMMA - 1)) - 1))

    # -- A0/A --
    @property
    def A_A0(self):
        return self._A_A0
    
    @A_A0.setter
    def A_A0(self, AA0: float, sup = True):
            
        sys1 = optimize(lambda M: AA0Mach(M, self.GAMMA))

        if sup:
            sys1.addGeq(1)
            ans = sys1.target(1 + 1e-5, AA0)
        
        else:
            sys1.addLeq(1)
            ans = sys1.target(1e-5, AA0)
            
        self.mach = ans
    
    
    # ----------------------------------
    # Constructors
    # ----------------------------------
    
    # --- temperature ---
    @staticmethod
    def from_temp(T0T: float | None = None, T0: float | None = None, T: float | None = None, GAMMA: Config | float = config):
        if isinstance(GAMMA, Config): GAMMA = GAMMA.GAMMA

        if T0T is None and (T0 is None or T is None):
            raise ValueError("Not enough inputs (Isen.from_temp)")
        
        elif T0T is None:
            T0T = T0 / T
        
        mach = np.sqrt(2 / (GAMMA - 1) * (T0T - 1))
        return Isen(mach, T = T, T0 = T0, GAMMA = GAMMA)
    
    # --- pressure ---
    @staticmethod
    def from_pres(P0P: float | None = None, P0: float | None = None, P: float | None = None, GAMMA: Config | float = config):
        if isinstance(GAMMA, Config): GAMMA = GAMMA.GAMMA
        
        if P0P is None and (P0 is None or P is None):
            raise ValueError("Not enough inputs (Isen.from_pres)")
        
        elif P0P is None:
            P0P = P0 / P
            
        mach = np.sqrt(2/(GAMMA - 1) * (np.pow(P0P, (GAMMA - 1)/(GAMMA)) - 1))
        
        return Isen(mach, P = P, P0 = P, GAMMA = GAMMA)
    
    # --- density ---
    @staticmethod
    def from_dens(r0_r: float | None = None, r0: float | None = None, r: float | None = None, GAMMA: Config | float = config):
        if isinstance(GAMMA, Config): GAMMA = GAMMA.GAMMA
        
        if r0r is None and (r0 is None or r is None):
            raise ValueError("Not enough inputs (Isen.from_dens)")
        
        elif r0r is None:
            r0r = r0 / r
            
        mach = np.sqrt(2/(GAMMA - 1) * (np.pow(r0r, (GAMMA - 1)/(GAMMA)) - 1))
        
        return Isen(mach, r = r, r0 = r0, GAMMA = GAMMA)
        
    # --- area ---
    @staticmethod
    def from_area(AA0: float | None = None, A0: float | None = None, A: float | None = None, GAMMA: Config | float = config, sup = True):
        if isinstance(GAMMA, Config): GAMMA = GAMMA.GAMMA
        
        if AA0 is None and (A0 is None or A is None):
            raise ValueError("Not enough inputs (Isen.from_area)")
        
        elif AA0 is None:
            AA0 = A / A0
        
        if (config.Q_(AA0).to('').m < 1):
            raise ValueError(f"Area Ratio is invalid (Isen.from_area) [{AA0:.3f} < 1]")
        
        
        
        sys = optimize(lambda M: AA0Mach(M, GAMMA))

        if sup:
            sys.addGeq(1)
            mach = sys.target(1 + 1e-5, AA0)
        
        else:
            sys.addLeq(1)
            mach = sys.target(1e-5, AA0)
        
        return Isen(mach, A = A, A0 = A0, GAMMA = GAMMA)
    
    # --- vel ---
    @staticmethod
    def from_vel(temp: float, vel: float, GAMMA: Config | float = config):
        if isinstance(GAMMA, Config): GAMMA = GAMMA.GAMMA
        
        mach = SpeedOfSound.mach(temp, vel)
        
        return Isen(mach, T = temp, GAMMA = GAMMA)
    
    @staticmethod
    def from_nu(nu2, GAMMA: Config | float = config, **kwargs):
        if isinstance(GAMMA, Config): GAMMA = GAMMA.GAMMA
        
        sys = optimize(lambda M: nuMach(M, GAMMA))
        
        sys.addGeq(1)
        mach = sys.target(1 + 1e-5, nu2)
        
        return Isen(mach, **kwargs)
    
    
    
    # ----------------------------------
    # Attribute Propagator
    # ----------------------------------
    
    def prop_up(rel: Relation):
        
        if rel.state1.T is not None:
            rel.state2.T = rel.state1.T * rel.T2_T1
            # config.safe_set(rel.state2, "T", rel.state1.T * rel.T2_T1)
       
        if rel.state1.P is not None:
            rel.state2.P = rel.state1.P * rel.P2_P1
            # config.safe_set(rel.state2, "P", rel.state1.P * rel.P2_P1)
           
        if rel.state1.A is not None:
            rel.state2.A = rel.state1.A * rel.A2_A1
            # config.safe_set(rel.state2, "A", rel.state1.A * rel.A2_A1)
            
        from .Relations.normal import Normal
        if isinstance(rel, Normal):
           
            if rel.state1.mach is not None:
                
                M2_num = np.pow(rel.state1.mach, 2) * (rel.GAMMA - 1) + 2
                M2_den = 2*rel.GAMMA*np.pow(rel.state1.mach, 2) - (rel.GAMMA - 1)
               
                mach2 = np.sqrt(M2_num / M2_den)
               
                print(mach2)
               
                config.safe_set(rel.state2, "mach", mach2)
      
    def prop_down(rel: Relation):
        
        if rel.state2.T is not None:
            rel.state1.T = rel.state2.T / rel.T2_T1
            # config.safe_set(rel.state1, "T", rel.state2.T / rel.T2_T1)
       
        if rel.state2.P is not None:
            rel.state1.P = rel.state2.P / rel.P2_P1
            # config.safe_set(rel.state1, "P", rel.state2.P / rel.P2_P1)
           
        if rel.state2.A is not None:
            rel.state1.A = rel.state2.A / rel.A2_A1
            # config.safe_set(rel.state1, "A", rel.state2.A / rel.A2_A1)
            
        from .Relations.normal import Normal
        if isinstance(rel, Normal):
           
            if rel.state2.mach is not None:
                num = np.sqrt((rel.GAMMA - 1)*np.pow(rel.state2.mach, 2) + 2)
                den = np.sqrt(rel.GAMMA * (2 * np.pow(rel.state2.mach, 2) - 1) + 1)
    
                M1 = num / den                
                config.safe_set(rel.state1, "mach", M1)

    
    def attr_pulse(self):
        self._propagate()
    
    def _propagate(self, source = None):
            
        if self._updating: return
        self._updating = True
        
        try:
            
            for child, rel in getattr(self, '_children', []):
                if child is not source:
                    Isen.prop_up(rel)
                    child._propagate(source = self)
            
            for parent, rel in getattr(self, '_parent', []):
                if parent is not source:
                    Isen.prop_down(rel)
                    parent._propagate(source = self)
            
        finally:
            self._updating = False
    
    # OLD METHOD
    # @staticmethod # 1 -> 2
    # def prop_up(rel: Relation, visited):
    #    
    #     if rel.state1.T is not None:
    #         config.safe_set(rel.state2, "T", rel.state1.T * rel.T2_T1)
    #    
    #     if rel.state1.P is not None:
    #         config.safe_set(rel.state2, "P", rel.state1.P * rel.P2_P1)
    #        
    #     if rel.state1.A is not None:
    #         config.safe_set(rel.state2, "A", rel.state1.A * rel.A2_A1)
    #        
    #     from .Relations.normal import Normal
    #     if isinstance(rel, Normal):
    #        
    #         if rel.state1.mach is not None:
    #            
    #             M2_num = np.pow(rel.state1.mach, 2) * (rel.GAMMA - 1) + 2
    #             M2_den = 2*rel.GAMMA*np.pow(rel.state1.mach, 2) - (rel.GAMMA - 1)
    #            
    #             mach2 = np.sqrt(M2_num / M2_den)
    #            
    #             config.safe_set(rel.state2, "mach", mach2)
    #        
    #        
    #        
    #
    # @staticmethod # 2 -> 1
    # def prop_down(rel: Relation, visited):
    #    
    #     if rel.state2.T is not None:
    #         config.safe_set(rel.state1, "T", rel.state2.T / rel.T2_T1)
    #    
    #     if rel.state2.P is not None:
    #         config.safe_set(rel.state1, "P", rel.state2.P / rel.P2_P1)
    #        
    #     if rel.state2.A is not None:
    #         config.safe_set(rel.state1, "A", rel.state2.A / rel.A2_A1)
    #            
    #     from .Relations.normal import Normal
    #     if isinstance(rel, Normal):
    #        
    #         if rel.state2.mach is not None:
    #             num = np.sqrt((rel.GAMMA - 1)*np.pow(rel.state2.mach, 2) + 2)
    #             den = np.sqrt(rel.GAMMA * (2 * np.pow(rel.state2.mach, 2) - 1) + 1)
    #
    #             M1 = num / den                
    #             config.safe_set(rel.state1, "mach", M1)
    #
    #
    # def attr_pulse(self):
    #     visited = set()
    #     active = set()
    #     self._propagate(visited, active)
    #     return visited
    #
    # def _propagate(self, visited, active):
    #   
    #     if self.UUID in visited:
    #         return
    #   
    #     if self.UUID in active:
    #         return
    #   
    #     active.add(self.UUID)        
    #   
    #     for parent, rel in self._parent: # propagate upstream
    #         Isen.prop_up(rel, visited)
    #         parent._propagate(visited, active)
    #   
    #     for child, rel in self._children: # propagate downstream
    #         Isen.prop_down(rel, visited)
    #         child._propagate(visited, active)
    #
    #     active.remove(self.UUID)
    #     visited.add(self.UUID)
    # 
    # 
    # 
    # """
    
    
    # ----------------------------------
    # Attribute Accessor / Modifier
    # ----------------------------------
            
    def __setattr__(self, name: str, value):
        
        _variables = (
            "mu", "nu", "vel",
            "Tstar_T", "Pstar_P", "rstar_r", "Astar_A", "mach_star"
        )
        
        _internal_var = (
            "_mach", "GAMMA", "UUID",
            "_parent", "_children", "_tracker", "_updating", "_generating",
            "_T0_T", "_P0_P", "_r0_r", "_A_A0",
        )
        
        updating = False
        
        # If name is an internal variable or is an external variable while generating
        if (name in _internal_var) or (name in _variables and getattr(self, "_generating", False)):
            super().__setattr__(name, value)
            return
        
        m = re.match("^((T|P|r|A)0?)", name)
        
        if m: # T(0) / P(0) / r(0) / A(0)
            var, i = m.groups()

            # if generating and already defined, then pass
            if getattr(self, "_generating", False) and (getattr(self, var, None) is not None): return
            
            super().__setattr__(name, value)
            
            if value is None: return
            
            try:
            
                match var:    
                    # case "T":
                    #     super().__setattr__("T0", value * self.T0_T)
                    #     super().__setattr__("Tstar", value * self.Tstar_T)
                    
                    case "T0":
                        self.T = value / self.T0_T
                        
                    # case "P":
                    #     super().__setattr__("P0", value * self.P0_P)
                    #     super().__setattr__("Pstar", value * self.Pstar_P)
                                                
                    case "P0":
                        self.P = value / self.P0_P

                    # case "r":
                    #     super().__setattr__("r0", value * self.r0_r)
                    #     super().__setattr__("rstar", value * self.rstar_r)
                    
                    case "r0":
                        self.r = value / self.r0_r
                        
                    case "A":
                        super().__setattr__("A0", value / self.A_A0)
                    
                    case "A0":
                        self.A = value * self.A_A0
                   
            finally:
                pass
        
        
        # if T is defined, calculate velocity
        if getattr(self, "T", None) is not None:
            super().__setattr__("vel", SpeedOfSound.vel(self.T, self.mach))
        
        validIdealGas = [i for i in ("T", "P", "r") if getattr(self, f"{i}", None) is None]
        
        ## If one is undefined, then fill the missing value
        if len(validIdealGas) <= 1:
            err = lambda E, T: abs(E - T) / E
            
            if len(validIdealGas) == 1:
                match validIdealGas[0]:
                    case "T":
                        self.T = IdealGas.temp(self.P, self.r)
                    case "P":
                        self.P = IdealGas.pres(self.T, self.r)
                    case "r":
                        self.r = IdealGas.dens(self.P, self.T)
                        
            isInvalid = err(IdealGas.dens(self.P, self.T).to_base_units().m, config.Q_(self.r).to_base_units().m) > 1e-5
            
            if isInvalid:
                if name == "r":
                    self.T = IdealGas.temp(self.P, self.r)
                else:
                    self.r = IdealGas.dens(self.P, self.T)
        
        self.attr_pulse()
        
        # If name is already defined
        if name in dir(self):
            super().__setattr__(name, value)
            return
        
        raise AttributeError(f"Invalid assignment, Isen has no variable [{name}]")
    
    
    def __getattr__(self, name):
        
        if name not in ["mach", "_mach"]: M = config.Q_(self.mach).to('').m
        
        match name:
            case "nu":
                if M >= 1:
                    nu = nuMach(self.mach)
                    return nu.to('deg')
            
            case "mu":
                if M >= 1:
                    mu = np.asin(1 / self.mach)
                    return config.Q_(mu).to('deg')
                
            case "A0_A":
                return 1 / self.A_A0                

        m = re.match(r"^(T|P|r)_\1(0)$", name)
        
        if m:
            var, i = m.groups()
            return 1 / getattr(self, f"{var}0_{var}")
        
        m = re.match(r"^(T|P|r)(0|star)$", name)
        
        if m:
            var, i = m.groups()
            
            if getattr(self, f'{var}', None) is None:
                return None
            
            if i == '0':
                return getattr(self, f"{var}") * getattr(self, f"{var}0_{var}")

            else:
                return getattr(self, f"{var}") * getattr(self, f"{var}star_{var}")
        
        
        raise AttributeError(f"Isen does not have this attribute [{name}]")
    
    
    # ----------------------------------
    # String
    # ----------------------------------
    def __str__(self):
        
        line1 = "="*4 + " Isentropic "+"="*4 + f"\nmach: {self.mach}" + "\n"
        
        line1b = "Mach angle (mu): " + (f"{self.mu}" if (self.mach.to('').m >= 1) else "[n/a]")
        line1c = "P-M angle (nu): " + (f"{self.nu}"  if (self.mach.to('').m >= 1) else "[n/a]")
        
        line1 += f"{line1b}\t{line1c}\n" + "-"*20 + "\n"
                
        line2 = f"T0/T: {self.T0_T}\tP0/P: {self.P0_P}\tr0/r: {self.r0_r}\n" + "-"*20 + "\n"
        
        line3 = f"T*/T: {self.Tstar_T}\tP*/P: {self.Pstar_P}\n"
        line4 = f"r*/r: {self.rstar_r}\tA/A*: {self.A_A0}\n" + "-"*20 + "\n"
        
        Pstr = f"P: {self.P}\t P0: {self.P0}\t Pstar: {self.Pstar}\n" if (getattr(self, "P", None) is not None) else ""
        Tstr = f"T: {self.T}\t T0: {self.T0}\t Tstar: {self.Tstar}\n" if (getattr(self, "T", None) is not None) else ""
        rstr = f"r: {self.r}\t r0: {self.r0}\t rstar: {self.rstar}\n" if (getattr(self, "r", None) is not None) else ""
        Astr = f"A: {self.A}\t A0: {self.A0}\n" if (getattr(self, "A", None) is not None) else ""
        
        return line1 + line2 + line3 + line4 + Pstr + Tstr + rstr + Astr + "="*20 + "\n"



def add_rel(state_parent: Isen, state_child: Isen, rel = None):
    state_parent._parent.append((state_child, rel))
    state_child._children.append((state_parent, rel))
