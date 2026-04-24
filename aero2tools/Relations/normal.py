from .relation import *

# ================================
# Normal Shock
# ================================
"""
IN:
    P2/P1
    r2/r1
    T2/T1
    P02/P01
    P1/P02
"""

class Normal(Relation):
    
    def __init__(self, state1: Isen):
        super().__init__()
        
        self.state1 = state1
        
        M2_num = np.pow(state1.mach, 2) * (self.GAMMA - 1) + 2
        M2_den = 2*self.GAMMA*np.pow(state1.mach, 2) - (self.GAMMA - 1)
        
        mach2 = np.sqrt(M2_num / M2_den)
        
        self.state2 = Isen(mach2)
        
        self.state1._children.append((self.state2, self))
        self.state2._parent.append((self.state1, self))
        
        self.state1.attr_pulse()
        
    def propagate(self, forward = True):
        changed = False
        
        s1 = self.state1
        s2 = self.state2        
        
        if forward and s1.mach is not None:
            g = self.GAMMA
            M1 = s1.mach
            
            M2 = np.sqrt((M1*M1*(g-1)+2)/(2*g*M1*M1-(g-1)))

            if not config.approx(s2.mach, M2):
                s2.mach = M2
                changed = True
                
        if not forward and s2.mach is not None:
            g = self.GAMMA
            M2 = s2.mach
            
            M1 = np.sqrt((M2*M2*(g-1)+2)/(2*g*M2*M2-(g-1)))

            if not config.approx(s1.mach, M1):
                s1.mach = M1
                changed = True
        
        changed = super().propagate(forward) or changed
        
        return changed
    
    # ------------------
    # constructors
    # ------------------
    
    @staticmethod
    def from_temp(T2_T1 = None, T1 = None, T2 = None):
        
        [T2_T1, T2, T1] = config.ratio_check(T2_T1, T2, T1, 'Normal.from_temp')
                
        state1 = Isen(2, T = T1)
        nrm = Normal(state1)
        nrm.T2_T1 = T2_T1
        
        return nrm
    
    @staticmethod
    def from_pres(P2_P1 = None, P1 = None, P2 = None):
        
        if (P2_P1 is None) and (P1 is None or P2 is None):
            raise ValueError("Not enough inputs (Normal.from_pres)")
        elif (P2_P1 is None):
            P2_P1 = P2 / P1
        
        state1 = Isen(2, P = P1)
        nrm = Normal(state1)
        nrm.P2_P1 = P2_P1
        
        return nrm
    
    @staticmethod
    def from_dens(r2_r1 = None, r1 = None, r2 = None):
        
        if (r2_r1 is None) and (r1 is None or r2 is None):
            raise ValueError("Not enough inputs (Normal.from_dens)")
        elif (r2_r1 is None):
            r2_r1 = r2 / r1
        
        state1 = Isen(2, r = r1)
        nrm = Normal(state1)
        nrm.r2_r1 = r2_r1
        
        return nrm
    
    
    
    
    # ------------------
    # properties
    # ------------------
    
    # --- density ---
    @property
    def r2_r1(self):
        
        num = (self.GAMMA + 1) * np.pow(self.state1.mach, 2)
        den = (self.GAMMA - 1)*np.pow(self.state1.mach, 2) + 2
        
        return num / den
    
    @r2_r1.setter
    def r2_r1(self, r2_r1):
        
        num = 2 * r2_r1
        den = (self.GAMMA + 1) - r2_r1 * (self.GAMMA - 1)
        
        M = np.sqrt(num / den)
        self.state1.mach = M
        
    # --- pressure ---
    @property
    def P2_P1(self):
        
        num = 2*self.GAMMA*np.pow(self.state1.mach, 2) - (self.GAMMA - 1)
        den = self.GAMMA + 1
        
        return num / den
    
    @P2_P1.setter
    def P2_P1(self, P2_P1):
        
        num = P2_P1 * (self.GAMMA + 1) + (self.GAMMA - 1)
        den = 2 * self.GAMMA
        
        M = np.sqrt(num / den)
        self.state1.mach = M
        
    # --- temperature ---
    @property
    def T2_T1(self):
        return self.P2_P1 / self.r2_r1
    
    @T2_T1.setter
    def T2_T1(self, T2_T1):
        if (T2_T1 <= 1): raise ValueError("T2/T1 must be greater than 1 for Normal Shocks")
        
        g = self.GAMMA
        x = T2_T1
        
        a = ( g**2 * x**2 + 2*g**2 * x + g**2 + 2*g * x**2 - 12*g * x + 2*g + x**2 + 2*x + 1 ) 
        
        # Full numerator and denominator of the Mach expression 
        num = (g + 1) * np.sqrt(a) + (g**2 + 1) * (x + 1) + 2*g*x - 6*g 
        den = 4 * (g - 1) * g 
        
        M = np.sqrt(num / den)
        self.state1.mach = M
        
    # --- stagnation pressure ---
    @property
    def P02_P01(self):
        return self.P2_P1 * np.pow(self.T2_T1, -(self.GAMMA)/(self.GAMMA - 1))
    
    # --- pitot pressure ---
    @property
    def P02_P1(self):
        M = self.state1.mach
        
        t1 = np.pow((self.GAMMA + 1)*M, 2)/(4*self.GAMMA*np.pow(M, 2) - 2*(self.GAMMA - 1))
        t2 = ((1 - self.GAMMA) + 2*self.GAMMA*np.pow(M,2)) / (self.GAMMA + 1)
        
        return np.pow(t1, self.GAMMA/(self.GAMMA - 1)) * t2
    
    @P02_P1.setter
    def P02_P1(self, P02_P1):
        
        num = 9*np.pow(self.GAMMA + 1, 2)
        den = 2*(17*self.GAMMA + 1)        
        
        K = np.pow(num/den, self.GAMMA / (self.GAMMA - 1))
        
        
        a_num = 2*self.GAMMA*(289*self.GAMMA + 19)
        a_den = 9*(self.GAMMA + 1)*(17*self.GAMMA+1)
        a1 = a_num / a_den * K
        
        b1 = (34*self.GAMMA)/(3*(self.GAMMA+1)) * K
        
        c1 = (17*self.GAMMA + 1)/(self.GAMMA + 1) * K
        
        a = a1 / 2
        b = b1 - 3*a1
        c = 9/2*a1 - 3*b1 + c1
        
        M = (-b + np.sqrt(np.pow(b, 2) - 4*a*(c - P02_P1)))/2/a
        self.state1.mach = M
    
    @property
    def P1_P02(self):
        return 1 / self.P02_P1
    
    # ----------------------------------
    # String
    # ----------------------------------
    def __str__(self):
        
        line1 = "="*4 + " Normal "+"="*4 + "\n"
        line2 = f"M1: {self.mach1}\tM2: {self.mach2}\n"
        line3 = f"P02/P01: {self.P02_P01}\tP1/P02: {self.P1_P02}\n"
        line4 = f"P2/P1: {self.P2_P1}\tr2/r1: {self.r2_r1}\tT2/T1: {self.T2_T1}\n"
        line5 = "="*20 + "\n"
        return line1 + line2 + line3 + line4 + line5