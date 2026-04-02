from .core import *
from .isentropic import *
from .shock import *

import numpy as np


class Airfoil:
    
    def __init__(self, state0, alpha, c1, c2, tau):
        
        self.state0 = state0
        
        self.c1 = c1 + np.finfo(np.float32).eps
        self.c2 = c2 + np.finfo(np.float32).eps
        self.tau = tau
        
        self.theta1 = np.atan(self.tau / self.c1 / 2)
        self.theta2 = np.atan(self.tau / self.c2 / 2)
        
        self.alpha = alpha
    
    def FlatPlate(state0, alpha = config.Q_(0, 'deg')):
        return Airfoil(state0, alpha, 1, 0, 0)
    
    def Kite(state0, theta1 = config.Q_(0, 'deg'), theta2 = config.Q_(0, 'deg'), alpha = config.Q_(0, 'deg')):
        
        if (theta1.m == 0 and theta2.m == 0):
            c2 = 0
        else:
            c2 = np.tan(theta1) / (np.tan(theta1) + np.tan(theta2))
        
        c1 = 1 - c2
        tau = 2*np.tan(theta1)*c1
        
        return Airfoil(state0, alpha, c1, c2, tau)
    
    @property
    def alpha(self):
        return self._alpha
    
    @alpha.setter
    def alpha(self, alpha):
        self._alpha = alpha
        
        self.top = IsenTracker(self.state0)
        self.top.addDeflection((self.theta1 - alpha))
        self.top.addDeflection(-(self.theta1 + self.theta2))
        
        self.bot = IsenTracker(self.state0)
        self.bot.addDeflection((self.theta1 + alpha))
        self.bot.addDeflection(-(self.theta1 + self.theta2))
    
    @property
    def comp(self):
        
        Cx = (self.P3_Pinf + self.P1_Pinf)*self.tau - (self.P4_Pinf + self.P2_Pinf)*self.tau
        Cx /= config.GAMMA * self.mach**2
                
        Cy = (self.P3_Pinf - self.P1_Pinf)*self.c1 + (self.P4_Pinf - self.P2_Pinf)*self.c2
        Cy /= config.GAMMA * self.mach**2 / 2
        
        return [Cx, Cy]
    
    @property
    def cd(self):
        Cx, Cy = self.comp
        return Cx * np.cos(self.alpha) + Cy * np.sin(self.alpha)
    
    @property
    def cl(self):
        Cx, Cy = self.comp
        return Cy * np.cos(self.alpha) - Cx * np.sin(self.alpha)
    
    
    def __getattr__(self, name):
        
        try:
            
            n_top = len(self.top.states) - 1
            n_bot = len(self.bot.states) - 1

            
            m = re.match(rf"^(mach|{'|'.join([f'{i}0' for i in IsenTracker._variables_A])}|{'|'.join(IsenTracker._variables_A)})([0-9]*)$", name)
            
            if m:
                var, i = m.groups()
                i = (int(i)) if (i != '') else 0
                
                if i <= 0:
                    return getattr(self.top.states[0], var)
                
                elif (i <= n_top):
                    var = f"{var}2"
                    return getattr(self.top.states[i], var)
                
                elif (i <= n_bot  + n_top):
                    i -= n_top
                    var = f"{var}2"
                    return getattr(self.bot.states[i], var)
    
    
            m = re.match(rf"^P([0-9]*)_Pinf$", name)
            
            if m:
                i = m.groups()
                i = (int(i[0])) if (i[0] != '') else 0
                
                if i == 0:
                    return config.Q_(1, '')
                
                elif (i <= n_top):
                    return getattr(self.top, f'P{i + 1}_P1')
                
                elif (i <= n_top + n_bot):               
                    i -= n_top
                    return getattr(self.bot, f'P{i + 1}_P1')
                
    
        except AttributeError:
            pass
        
        raise AttributeError(f"Kite has no attribute ({name})")