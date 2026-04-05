from .core import *






class optimize:
    
    h = 1e-5
    epsilon = 1e-10
    
    def __init__(self, func):
        
        
        self.EQ = []
        self.INEQ = []
        
        self.func = func
        self.rho = 1
        
    
    def funcp(self, x, p):
        
        return self.func(x) + \
            sum([(1) if (g[0](x) > 0) else (0) for g in self.INEQ]) + \
            p * sum([max(g[0](x), 0)**2 for g in self.INEQ])
            
    def addLeq(self, val):
        self.addIneq(val, False)
    
    def addGeq(self, val):
        self.addIneq(val, True)
        
    def addIneq(self, val, flipped = False):
        '''flipped: True [g(x) >= 0], False [g(x) <= 0]'''
        if not flipped:
            def tempfunc(x):
                return x - val
        
        else:
            def tempfunc(x):
                return val - x
            
        self.INEQ.append((tempfunc, val, flipped))         
    
        
    def optimize(self, init):
        
        p = 1
        xN = init
        
        for j in range(190):
            d1 = diff(lambda x: self.funcp(x, p), xN)
            d2 = diff2(lambda x: self.funcp(x, p), xN)
            
            dx = d1 / d2
            xN -= dx
            
            if np.abs(dx) < optimize.epsilon: break
            p *= 1.2
        
        return xN
    
    def target(self, init, targ):
        
        p = 1
        xN = init
        
        for j in range(100):
            f = self.funcp(xN, p)
            d = diff(lambda x: self.funcp(x, p), xN)
            
            dx = (f - targ) / d
            xN -= dx
            
            if np.abs(dx) < optimize.epsilon: break
            p *= 1.2
        
        return xN
    
    
    
def diff(func, x):
    
    f1 = func(x)
    f2 = func(x + optimize.h)
    
    return (f2 - f1) / optimize.h

def diff2(func, x):
    
    f1 = diff(func, x)
    f2 = optimize.diff(func, x + optimize.h)
    
    return (f2 - f1) / optimize.h