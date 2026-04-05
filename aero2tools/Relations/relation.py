from __future__ import annotations
from ..core import *
from ..optimize import *
from ..state import *

import re
import numpy as np


class Relation:
    
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
    def P1_P02(self):
        return self.state1.P_P0 * self.P02_P01
    
    @property
    def entropy(self):
        return config.CP * np.log(self.T2_T1) - config.R * np.log(self.P2_P1)
    
    