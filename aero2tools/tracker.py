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
        
        
        
        
        
        
        
















