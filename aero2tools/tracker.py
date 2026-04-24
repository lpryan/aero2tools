from .core import *
from .optimize import *
from .state import *

from .Relations.interIsen import InterIsen
from .Relations.normal import Normal
from .Relations.oblique import Oblique

from abc import ABC, abstractmethod
import numpy as np
import re

# ===========================
# Isentropic State Tracker
# ===========================

_state_var_ = (
    "mach", "nu", "mu",
)

_state_var_ratio_ = (
    "T", "P", "r", "A"
)

_relation_var_ = (
    "fwd", "rwd",
)


class Tracker:
    
    def __init__(self, state0: Isen):
        
        self.states = [state0]
        self.relations = []
        
    def addRelation(self, func, theta = None, **kwargs):
        
        inter = None
        statePrev = self.states[-1]
        
        if func is Normal:
            inter = Normal(statePrev)
        
        elif func is Oblique:       
            
            # check if theta or beta are defined
            
            if theta is not None:
                inter = Oblique.from_theta(statePrev, theta, False)
    
        elif func is InterIsen:
            
            if theta is not None:
                inter = InterIsen.from_theta1(statePrev, theta)
            
        else:
            raise ValueError("Invalid Relation type (addRelation)")
    
        if inter is None:
            raise ValueError("Not enough inputs (addRelation)")
        
        self.states.append(inter.state2)
        self.relations.append(inter)
        
    def __getattr__(self, name):
        
        if name in ["states", "relations"]:
            return object.__getattr__(self, name)
        
        m = re.match(rf"^({'|'.join([f'{i}0' for i in _state_var_ratio_])}|{'|'.join(_state_var_ratio_)}|{'|'.join(_state_var_)})([0-9]*)$", name)
        
        if m:
            var, i = m.groups()
            i = (int(i)) if (i != '') else 0
            
            if i != 0: i -= 1
            
            return getattr(self.states[i], var)
        
        m = re.match(rf"^({'|'.join(_relation_var_)})([0-9]*)$", name)
        
        if m:
            var, i = m.groups()
            i = int(i) if (i != '') else 0
            
            if i != 0: i -= 1
            
            return getattr(self.relations[i], var)
        
        
        
        
        
        
        
        
        m = re.match(rf"^(P0|{'|'.join(_state_var_ratio_)})(\d+)_\1(\d+|)$", name)
        
        if m:
            var, i, j = m.groups()
            
            i = int(i) if (i != '') else 1
            j = int(j) if (j != '') else 1
            
            if (i == 0):
                if j >= len(self.states):
                    raise IndexError(f"Invalid assumption of number of states [{j} > {len(self.states)}]")
                                
                return getattr(self.states[j - 1], f"{var}0_{var}")
            
            elif (i - j == 1):
                if j >= len(self.states):
                    raise IndexError(f"Invalid assumption of number of states [{i} > {len(self.states)}]")
                
                return getattr(self.relations[j - 1], f"{var}2_{var}1")
            
            elif (j > i):
                return 1 / self.__getattr__(f"{var}{j}_{var}{i}")
            
            else:
                prod = 1
                
                for n in range(j, i):
                    prod *= self.__getattr__(f"{var}{n+1}_{var}{n}")
                    
                return prod
        
    
        raise AttributeError(f"Tracker does not have this attribute [{name}]")
    
    
    
    
    def __setattr__(self, name, value):
        
        if name in ["states", "relations"]:
            object.__setattr__(self, name, value)
            return
        
        
        m = re.match(rf"^({'|'.join([f'{i}0' for i in _state_var_ratio_])}|{'|'.join(_state_var_ratio_)}|{'|'.join(_state_var_)})([0-9]*)$", name)
        
        if m:
            var, i = m.groups()
            i = (int(i)) if (i != '') else 0
            
            if i != 0: i -= 1
            
            return setattr(self.states[i], var, value)
        
        raise AttributeError(f"Invalid assignment, Tracker has no variable [{name}]")
        
        
















