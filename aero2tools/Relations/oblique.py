from .relation import *

import numpy as np

# ================================
# Oblique Shock
# ================================
"""
IN:
    M1 +
    
    turn angle (weak) or
    turn angle (strong) or
    wave angle or
    mach1n

OUT:
    M2, turn angle, wave angle
    P2/P1, r2/r1, T2/T1,
    P02/P01, mach1n, mach2n
"""

def thetaMachBeta(mach1, beta):
    
    num = np.pow(mach1*np.sin(beta), 2) - 1
    den = mach1*mach1*(config.GAMMA + np.cos(2*beta)) + 2
    
    tan_theta = 2/np.tan(beta) * num / den
    return np.atan(tan_theta)

def BetaMax(mach):
    
    sys = optimize(lambda B: thetaMachBeta(mach, B))
    sys.addLeq(np.pi / 2)
    sys.addGeq(0)    
    
    mach = config.Q_(mach).to('').m
    beta_max = sys.optimize(np.pi / 4)
    return config.Q_(beta_max, 'rad')





class Oblique(Relation):
    
    def __init__(self, state1: Isen, beta):
        super().__init__()
        
        self.state1 = state1
        self.beta = beta
        
        self.state1._children.append((self.state2, self))
        self.state2._parent.append((self.state1, self))
        
        self.state1.attr_pulse()
    
    
    def propagate(self, forward = True):
        changed = False
        
        s1 = self.state1
        s2 = self.state2        
        
        if forward and s1.mach is not None:
            
            if not config.approx(s2.mach, self.mach2):
                s2.mach = self.mach2
                changed = True
                
        if not forward and s2.mach is not None:
            g = self.GAMMA
            M2 = s2.mach
            B = self.beta
            
            pass
            # M1n = np.sqrt((M2n*M2n*(g-1)+2)/(2*g*M2n*M2n-(g-1)))
            # M1 = M1n / np.sin(B)

            # if s2.mach != M1:
            #     s1.mach = M1
            #     changed = True
        
        changed = super().propagate(forward) or changed
                
        return changed
    
    # ------------------
    # constructors
    # ------------------
    
    @staticmethod # strong (M2 < 1) | weak (M2 > 1)
    def from_theta(state1: Isen | float, theta, strong = True):
        Q_ = config.Q_
        
        if isinstance(state1, float):
            state1 = Isen(state1)
            
        mach1 = Q_(state1.mach).to('').m
        
        if (mach1 < 1):
            raise ValueError(f"mach is subsonic ({mach1} < 1)")
        
        h = 1e-20
        
        mu = state1.mu.to('rad').m
        theta_target = Q_(theta).to('rad').m
        
        beta1 = (np.pi/2 - h) if strong else (mu + h)
        
        beta_max = BetaMax(mach = mach1)
        theta_max = thetaMachBeta(mach1, beta_max)
        
        if theta_target > theta_max:
            raise ValueError(f"theta exceeds maximum ({theta:~.4g} > {theta_max:~.4g})\n[M1: {mach1}, βmax: {beta_max}]")
        
        elif theta_target < 0:
            raise ValueError(f"theta is less than zero ({theta:~4g} < 0)")
    
        sys = optimize(lambda B: thetaMachBeta(mach1, B))
        sys.addLeq(np.pi / 2)
        sys.addGeq(0)
        
        beta_opt = sys.target(beta1, theta_target)
        theta_opt = thetaMachBeta(mach1, beta_opt)
        
        if not config.approx(theta_opt, theta_target):
            raise TimeoutError("Could not converge theta")
        
        inter = Oblique(state1, beta_opt)
        return inter
    
    @staticmethod
    def from_thetaBeta(theta, beta):
        Q_ = config.Q_
        
        THETA = Q_(theta).to('rad').m
        BETA = Q_(beta).to('rad').m
                
        sys = optimize(lambda M: thetaMachBeta(M, BETA))
                
        mach = sys.target(10, THETA)
        
        print(mach)
        
        state1 = Isen(mach)        
        
        inter = Oblique(state1, BETA)
        return inter      
        
    
    
    
    
    
    # ------------------
    # properties
    # ------------------
    
    # --- beta ---
    @property
    def beta(self):
        return config.Q_(getattr(self, "_beta", None)).to('deg')
    
    @beta.setter
    def beta(self, B):
        angle = config.Q_(B).to('rad').m
    
        if (angle > np.pi / 2) or (angle < self.state1.mu.to('rad').m):
            raise ValueError("Shock is detached")
        
        setattr(self, "_beta", B)
        
        if getattr(self, 'state2', None) is None:
            self.state2 = Isen(self.mach2)
        else:
            self.state2.mach = self.mach2
            
    @property
    def beta_max(self):
        return BetaMax(self.mach)
    
    @property
    def theta_max(self):
        return thetaMachBeta(self.mach1, self.beta_max)
    
    # --- mach1n ---
    @property
    def mach1n(self):
        return self.state1.mach * np.sin(self.beta)
    
    @mach1n.setter
    def mach1n(self, M1n):
        B = np.asin(M1n / self.state1.mach)
        self.beta = B
    
    # --- mach1t ---
    @property
    def mach1t(self):
        return self.state1.mach * np.cos(self.beta)
    
    ## === Normal Components ===
    
    @property
    def T2_T1(self):
        return self.P2_P1 / self.r2_r1
    
    @property
    def P2_P1(self):
        num = 2*self.GAMMA*np.pow(self.mach1n, 2) - (self.GAMMA - 1)
        den = self.GAMMA + 1
        return num / den
    
    @property
    def r2_r1(self):
        num = (self.GAMMA + 1) * np.pow(self.mach1n, 2)
        den = (self.GAMMA - 1)*np.pow(self.mach1n, 2) + 2
        return num / den    
        
    @property
    def P02_P01(self):
        return self.P2_P1 * np.pow(self.T2_T1, -(self.GAMMA)/(self.GAMMA - 1))
    
    # --- mach2n ---
    @property
    def mach2n(self):
        
        M2_num = np.pow(self.mach1n, 2) * (self.GAMMA - 1) + 2
        M2_den = 2*self.GAMMA*np.pow(self.mach1n, 2) - (self.GAMMA - 1)
        
        mach2n = np.sqrt(M2_num / M2_den)
    
        return mach2n
    
    @property
    def mach2t(self):
        return self.mach1t * np.sqrt(1 / self.T2_T1)
    
    @property
    def mach2(self):
        return np.sqrt(self.mach2n**2 + self.mach2t**2)
    
    @property
    def theta(self):
        return config.Q_(thetaMachBeta(self.mach1, self.beta)).to('deg')
    
    
    
    
    
    
    
    def __setattr__(self, name, value):
        
        _internal_vars = ("_beta")
        
        if name in _internal_vars:
            object.__setattr__(self, name, value)
            return
        
        return super().__setattr__(name, value)
    
    # ----------------------------------
    # String
    # ----------------------------------
    def __str__(self):
        
        line1 = "="*4 + " Oblique "+"="*4 + "\n"
        line2 = f"M1: {self.mach1}\tM2: {self.mach2}\n"
        line3 = f"Turn ang. [θ]: {self.theta}\tWave ang. [β]: {self.beta}\n"  
        line4 = f"P2/P1: {self.P2_P1}\tr2/r1: {self.r2_r1}\tT2/T1: {self.T2_T1}\n"
        line5 = f"P02/P01: {self.P02_P01}\tM1n: {self.mach1n}\tM2n: {self.mach2n}\n"
        line6 = "="*20
        return line1 + line2 + line3 + line4 + line5 + line6
    