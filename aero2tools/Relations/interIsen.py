from .relation import *

# ================================
# Inter-Isentropic Flow Relations
# ================================
"""
IN:
    T2/T1 (1 or 2)
    P2/P1 (1 or 2)
    r1/r2 (1 or 2)
    A1/A2 (1 or 2)
    Δθ (1 or 2)
    (1 and 2)

P01 = P02 | T01 = T02 | r01 = r02 | A01 = A02
"""

class InterIsen(Relation):
    
    def __init__(self, state1: Isen, state2: Isen):
        super().__init__()
        
        self.state1 = state1
        self.state2 = state2
        
        self.state1._children.append((self.state2, self))
        self.state2._parent.append((self.state1, self))
        
        self.state1.attr_pulse()
        
    def propagate(self, forward = True):
        changed = False
        
        s1 = self.state1
        s2 = self.state2
        
        if forward and s1.A is not None:
            new = s1.A * self.A2_A1
            if not config.approx(s2.A, new):
                s2.A = new
                changed = True
                
        elif not forward and s2.A is not None:
            new = s2.A / self.A2_A1
            if not config.approx(s1.A, new):
                s1.A = new
                changed = True
        
        changed = super().propagate(forward) or changed

        return changed
    
        
    # ------------------
    # constructors
    # ------------------
    
    # --- Temperature ---
    @staticmethod
    def from_temp1(state1: Isen, T2_T1 = None, T1 = None, T2 = None):
        
        if T1 is not None: 
            state1.T = T1
        
        T1 = getattr(state1, 'T', T1)
        T0 = getattr(state1, 'T0', None)
        
        if (T2_T1 is None) and (T1 is None or T2 is None):
            raise ValueError("Not enough inputs (InterIsen.from_temp1)")
        
        elif not (T1 is None or T2 is None):
            T2_T1 = T2 / T1
            
        T0_T2 = state1.T0_T / T2_T1
        
        state2 = Isen.from_temp(T0_T2, T = T2)
        state2.T0 = T0
        
        return InterIsen(state1, state2)
    
    @staticmethod
    def from_temp2(state2: Isen, T2_T1 = None, T1 = None, T2 = None):
        
        if T2 is not None:
            state2.T = T2
        
        T2 = getattr(state2, 'T', T2)
        T0 = getattr(state2, 'T0', None)
        
        if (T2_T1 is None) and (T1 is None or T2 is None):
            raise ValueError("Not enough inputs (InterIsen.from_temp2)")
        
        elif not (T1 is None or T2 is None):
            T2_T1 = T2 / T1
            
        T0_T1 = state2.T0_T * T2_T1
        
        state1 = Isen.from_temp(T0_T1, T = T1)
        state1.T0 = T0
        
        return InterIsen(state1, state2)
    
    # --- Pressure ---
    @staticmethod
    def from_pres1(state1: Isen, P2_P1 = None, P1 = None, P2 = None):
        
        if P1 is not None: 
            state1.P = P1
        
        P1 = getattr(state1, 'P', P1)
        P0 = getattr(state1, 'P0', None)
        
        if (P2_P1 is None) and (P1 is None or P2 is None):
            raise ValueError("Not enough inputs (InterIsen.from_pres1)")
        
        elif not (P1 is None or P2 is None):
            P2_P1 = P2 / P1
            
        P0_P2 = state1.P0_P / P2_P1
        
        state2 = Isen.from_pres(P0_P2, P = P2)
        state2.P0 = P0
        
        return InterIsen(state1, state2)
    
    @staticmethod
    def from_pres2(state2: Isen, P2_P1 = None, P1 = None, P2 = None):
        
        if P2 is not None:
            state2.P = P2
        
        P2 = getattr(state2, 'P', P2)
        P0 = getattr(state2, 'P0', None)
        
        if (P2_P1 is None) and (P1 is None or P2 is None):
            raise ValueError("Not enough inputs (InterIsen.from_pres2)")
        
        elif not (P1 is None or P2 is None):
            P2_P1 = P2 / P1
            
        P0_P1 = state2.P0_P * P2_P1
        
        state1 = Isen.from_pres(P0_P1, P = P1)
        state1.P0 = P0
        
        return InterIsen(state1, state2)
    
    # --- Density ---
    @staticmethod
    def from_dens1(state1: Isen, r2_r1 = None, r1 = None, r2 = None):
        
        if r1 is not None: 
            state1.r = r1
        
        r1 = getattr(state1, 'r', r1)
        r0 = getattr(state1, 'r0', None)
        
        if (r2_r1 is None) and (r1 is None or r2 is None):
            raise ValueError("Not enough inputs (InterIsen.from_dens1)")
        
        elif not (r1 is None or r2 is None):
            r2_r1 = r2 / r1
            
        r0_r2 = state1.r0_r / r2_r1
        
        state2 = Isen.from_dens(r0_r2, r = r2)
        state2.r0 = r0
        
        return InterIsen(state1, state2)
    
    @staticmethod
    def from_dens2(state2: Isen, r2_r1 = None, r1 = None, r2 = None):
        
        if r2 is not None:
            state2.r = r2
        
        r2 = getattr(state2, 'r', r2)
        r0 = getattr(state2, 'r0', None)
        
        if (r2_r1 is None) and (r1 is None or r2 is None):
            raise ValueError("Not enough inputs (InterIsen.from_dens2)")
        
        elif not (r1 is None or r2 is None):
            r2_r1 = r2 / r1
            
        r0_r1 = state2.r0_r * r2_r1
        
        state1 = Isen.from_dens(r0_r1, r = r1)
        state1.r0 = r0
        
        return InterIsen(state1, state2)
    
    # --- Area ---
    @staticmethod
    def from_area1(state1: Isen, A2_A1 = None, A1 = None, A2 = None):
        
        if A1 is not None: 
            state1.A = A1
        
        A1 = getattr(state1, 'A', A1)
        A0 = getattr(state1, 'A0', None)
        
        if (A2_A1 is None) and (A1 is None or A2 is None):
            raise ValueError("Not enough inputs (InterIsen.from_area1)")
        
        elif not (A1 is None or A2 is None):
            A2_A1 = A2 / A1
            
        A2_A0 = A2_A1 / state1.A0_A
        
        state2 = Isen.from_area(A2_A0, A = A2)
        state2.A0 = A0
        
        return InterIsen(state1, state2)
    
    @staticmethod
    def from_area2(state2: Isen, A2_A1 = None, A1 = None, A2 = None):
        
        if A2 is not None:
            state2.A = A2
        
        A2 = getattr(state2, 'A', A2)
        A0 = getattr(state2, 'A0', None)
        
        if (A2_A1 is None) and (A1 is None or A2 is None):
            raise ValueError("Not enough inputs (InterIsen.from_area2)")
        
        elif not (A1 is None or A2 is None):
            A2_A1 = A2 / A1
            
        A1_A0 = A2_A1 * state2.A_A0
        
        state1 = Isen.from_area(A1_A0, A = A1)
        state1.A0 = A0
        
        return InterIsen(state1, state2)
    
    # --- theta (expansion) ---
    @staticmethod
    def from_theta1(state1: Isen, theta):
        
        nu1 = state1.nu
        nu2 = nu1 + theta
        
        state2 = Isen.from_nu(nu2)
        
        return InterIsen(state1, state2)
    
    @staticmethod
    def from_theta2(state2: Isen, theta):
        
        nu2 = state2.nu
        nu1 = nu2 - theta
        
        state1 = Isen.from_nu(nu1)
        
        return InterIsen(state1, state2)  
    
    
    # ------------------
    # properties
    # ------------------
    
    # --- pressure ---
    @property
    def P2_P1(self):
        return self.state1.P0_P / self.state2.P0_P
    
    @P2_P1.setter
    def P2_P1(self, P2P1):
        config.safe_set(self.state2, "P0_P", self.state1.P0_P / P2P1)
    
    # --- temperature ---
    @property
    def T2_T1(self):
        return self.state1.T0_T / self.state2.T0_T
    
    @T2_T1.setter
    def T2_T1(self, T2T1):
        config.safe_set(self.state2, "T0_T", self.state1.T0_T / T2T1)
        
    # --- density ---
    @property
    def r2_r1(self):
        return self.state1.r0_r / self.state2.r0_r
    
    @r2_r1.setter
    def r2_r1(self, r2r1):
        config.safe_set(self.state2, "r0_r", self.state1.r0_r / r2r1)
        
    # --- ---    
    @property
    def P02_P01(self):
        return 1
    
    @property
    def P1_P02(self):
        return self.state1.P_P0 * self.P02_P01
    
    # --- area ---
    @property
    def A2_A1(self):
        return self.state2.A_A0 / self.state1.A_A0
    
    def __getattr__(self, name):
        
        try:
            return super().__getattr__(name)
        
        except AttributeError:
            
            match name:
                
                case "theta":
                    return self.state2.nu - self.state1.nu
                
                case "fwd":
                    return self.state1.mu
                
                case "rwd":
                    return self.state2.mu - self.theta
            
        raise AttributeError(f"InterIsen does not have this attribute [{name}]")
        
        
        