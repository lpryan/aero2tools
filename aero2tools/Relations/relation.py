from ..core import *
from ..optimize import *
from ..state import *

import re
import numpy as np


class Relation:
    
    def __init__(self, GAMMA: Config | float = config):
        if isinstance(GAMMA, Config): GAMMA = GAMMA.GAMMA
        self.GAMMA = GAMMA
    
    def propagate(self, forward = True):
                
        changed = False
        
        s1 = self.state1
        s2 = self.state2
        
        
        if forward:
            
            if s1.T is not None:
                new = s1.T * self.T2_T1
                if not config.approx(s2.T, new):
                    s2.T = new
                    changed = True
            
            if s1.P is not None:
                new = s1.P * self.P2_P1
                if not config.approx(s2.P, new):
                    s2.P = new
                    changed = True
                    
            if s1.A is not None:
                new = s1.A * self.A2_A1
                if not config.approx(s2.A, new):
                    s2.A = new
                    changed = True
                        
        else: # not forward (backward)
            
            if s2.T is not None:
                new = s2.T / self.T2_T1
                if not config.approx(s1.T, new):
                    s1.T = new
                    changed = True
            
            if s2.P is not None:
                new = s2.P / self.P2_P1
                if not config.approx(s1.P, new):
                    s1.P = new
                    changed = True
                    
            if s2.A is not None:
                new = s2.A / self.A2_A1
                if not config.approx(s1.A, new):
                    s1.A = new
                    changed = True
        
        return changed
    
    
    @property
    def T2_T1(self):
        return None
    
    @property
    def P2_P1(self):
        return None
    
    @property
    def r2_r1(self):
        return None
        
    @property
    def P02_P01(self):
        return None
    
    @property
    def A2_A1(self):
        return self.state2.A_A0 / self.state1.A_A0
    
    
    @property
    def entropy(self):
        return config.CP * np.log(self.T2_T1) - config.R * np.log(self.P2_P1)
        
    def __getattr__(self, name):
        
        # get ratios, ex: P0/P1, A0/A1
        m = re.match(rf"^(P|T|r|A)(0|star)_\1(|1|2)$", name)
        
        if m:
            var, i, j = m.groups()
            if j == '': j = 1
            
            if i == 'star':
                match var:
                    case "P": 
                        return getattr(self, f"state{j}").Pstar_P
                    case "T": 
                        return getattr(self, f"state{j}").Tstar_T
                    case "r": 
                        return getattr(self, f"state{j}").rstar_r
                    
                    case "A":
                        return getattr(self, f"state{j}").Astar_A
            
            elif i == '0':
                match var:
                    case "P": 
                        return getattr(self, f"state{j}").P0_P
                    case "T": 
                        return getattr(self, f"state{j}").T0_T
                    case "r": 
                        return getattr(self, f"state{j}").r0_r
                    
                    case "A":
                        return getattr(self, f"state{j}").A0_A
        
        
        # get inverse ratios, ex: T1/T0, A2/A0
        m = re.match(rf"^(P|T|r|A)(1|2)_\1(0|star)$", name)
        
        if m:
            var, i, j = m.groups()
            return 1 / self.__getattr__(rf"{var}{j}_{var}{i}")
                
        
        # retrive state attributes, ex: P2, T1, mach1
        m = re.match(rf"^((P|T|r)(star|0)?|mach|vel|mu|nu)(|1|2)$", name)
        
        if m:
            var, i, j, k = m.groups()
            if k == '': k = 1
            return getattr(getattr(self, f"state{k}"), f"{var}", None)
        
        raise AttributeError(f"Relation does not have this attribute [{name}]")
        
        
    def __setattr__(self, name, value):
        
        _internal_vars = ("GAMMA", "state1", "state2", "_beta")
                
        if name in _internal_vars or is_property(self, name):
            super().__setattr__(name, value)
            return
        
        raise AttributeError(f"Invalid assignment, Relation has no variable [{name}]")
    