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
    "T", "P", "r"
)

_relation_var_ = (
    "T2_T1", "P2_P1", "r2_r1", "P02_P01", 
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
        
        m = re.match(rf"^({'|'.join([f'{i}0' for i in _state_var_ratio_])}|{'|'.join(_state_var_ratio_)}|{'|'.join(_state_var_)})([0-9]*)$", name)
        
        if m:
            var, i = m.groups()
            i = (int(i)) if (i != '') else 0
            
            if i != 0: i -= 1
            
            return getattr(self.states[i], var)
            
    
    
    def __setattr__(self, name, value):
        
        if name in ["states", "relations"]:
            object.__setattr__(self, name, value)
        
        
        m = re.match(rf"^({'|'.join([f'{i}0' for i in _state_var_ratio_])}|{'|'.join(_state_var_ratio_)}|{'|'.join(_state_var_)})([0-9]*)$", name)
        
        if m:
            var, i = m.groups()
            i = (int(i)) if (i != '') else 0
            
            if i != 0: i -= 1
            
            return setattr(self.states[i], var, value)
        
        
        
















