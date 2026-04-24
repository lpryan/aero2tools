from .relation import *

# ================================
# Inter-Isentropic Flow Relations
# ================================
"""
IN:
    Mach
    ToTo* (sub / sup)
    T/T* (above / below Tmax)
    P/P*
    PoPo* (sub / sup)
    U/U*
    (s* - s)/R (sub /  sup)

OUT:
    M
    To/To*
    T/T*
    P/P*
    Po/Po*
    U/U*
    (s* - s)/R
"""

class Rayleigh(Relation):
    
    def __init__(self, state1, state2):
        super().__init__()
        
        self.state1 = state1
        self.state2 = state2
        
        self.state1._children.append((self.state2, self))
        self.state2._parent.append((self.state1, self))
        
        self.state1.attr_pulse()
    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        