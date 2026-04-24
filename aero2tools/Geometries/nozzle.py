from ..tracker import *
from enum import Enum

"""
# Ae / At(hroat) -> determine the flow
# Pa / P0 -> determine sub/supersonic

# Subsonic -> highest mach is smallest area
# Supersonic -> lowest mach is smallest area

# choked if P0/Pa is greater than P/P* (@ mach = 1) 
# if choked flow -> Ae = A*
# otherwise A/A* = At/Ae * Ae/A*

#
# Enter Ae/At: ...
#
#
#
#
#
#

# Ae/At + Pa/Pt
   

# class Nozzle:
    
#     def __init__(self, AeAt, inter13, inter32, Pa = None, Pt0 = None):
        
#         self.AeAt = AeAt
        
#         self.interIsen = inter13
#         self.interNormal = inter32
        
#         self.crit1 = inter13.state1
#         self.crit3 = inter13.state2
        
#         self.crit2 = inter32.state2
        
#         self.Pa = Pa
#         self.Pt0 = Pt0
        
#         if not (Pa is None or Pt0 is None):
#             self.PbPt0 = Pa / Pt0
#         else: self.PbPt0 = None
        
        
#         self.Tt = None
    
    
#     @property
#     def type(self):
        
#         # if (self.Pa is None or self.Pt0 is None):
#         #     return NozzleShockPosition.subsonic
        
#         PbPt0 = self.PbPt0
                
#         if (PbPt0 > self.P1P0):
#             return NozzleShockPosition.subsonic
        
#         elif (self.P2P0 < PbPt0 < self.P1P0):
#             return NozzleShockPosition.shock_in_nozzle
        
#         elif (self.P3P0 < PbPt0 < self.P2P0):
#             return NozzleShockPosition.over_expanded
        
#         elif (PbPt0 < self.P3P0):
#             return NozzleShockPosition.under_expanded
    
#     @property
#     def Me(self):
        
#         # if (self.Pa is None or self.Pt0 is None):
#         #     return None
        
#         g = config.GAMMA
#         gm = g - 1
#         gp = g + 1
        
#         PbPt0 = self.PbPt0
#         AeAt = self.AeAt
        
#         mach1 = np.sqrt(2/gm * (np.pow(PbPt0, -gm/g) - 1))
        
#         mach2_sq = -1/gm + np.sqrt(1/gm/gm + 2/gm * np.pow(2/gp, gp/gm) / np.pow(AeAt * PbPt0, 2))
#         mach2 = np.sqrt(mach2_sq)
        
#         mach3 = self.crit3.mach
#         mach4 = self.crit3.mach
        
#         shock_type = self.type
        
#         if shock_type == NozzleShockPosition.subsonic:
#             return mach1
        
#         if shock_type == NozzleShockPosition.shock_in_nozzle:
#             return mach2
        
#         if shock_type == NozzleShockPosition.over_expanded:
#             return mach3
        
#         if shock_type == NozzleShockPosition.under_expanded:
#             return mach4
    
    
#     @property
#     def mdot(self):
        
#         g = config.GAMMA
#         P0 = None
#         T0 = None
#         Astar = None
                
#         match self.type:
#             case NozzleShockPosition.subsonic:
#                 P0    = self.crit1.P0
#                 T0    = self.crit1.T0
#                 Astar = self.crit1.Astar
            
#             case _:
#                 raise TypeError(f"Unsupported shock type ({self.type})")
        
#         if (P0 is None or T0 is None):
            
#             if (T0 is None):
#                 T0 = self.Tt
            
#             if (P0 is None):
#                 P0 = self.Pt0            
            
#             if (P0 is None or T0 is None):
#                 raise ValueError(f"invalid throat value")
        
#         mdot = np.sqrt(g / config.R / T0) * (5/6)**3 * P0 * Astar
        
#         return mdot
    
#     @property
#     def P1P0(self):
#         return self.crit1.P_P0
    
#     @property
#     def P2P0(self):
#         return self.crit3.P_P0 * self.interNormal.P2_P1
    
#     @property
#     def P3P0(self):
#         return self.crit3.P_P0
    
    
#     # =========================
#     # Constructor
#     # =========================
#     def Isen(AeAt = None, Ae = None, At = None, Pa = None, P0 = None):
        
#         if not (Ae is None or At is None):
#             AeAt = Ae / At
        
#         crit1 = Isen.from_area(AeAt, A0=At, A=Ae, sup = False)
#         crit3 = Isen.from_area(AeAt, A0=At, A=Ae, sup = True)
        
#         crit1.P = Pa
        
#         inter1 = InterIsen(crit1, crit3)
#         inter2 = Normal(crit3)
                
#         return Nozzle(AeAt, inter1, inter2, Pa = Pa, Pt0 = P0)
    
    
#     # P0 = 1 [atm]; Pe = 0.314 [atm]
#     # P0 = 5 [atm], At = 4.2 [in2], AeAt = 2.193
#     # AeAt = 1.616; Pe = 0.94 [atm]; P0 = 1 [atm]
#     # AeAt = 1.616; Pe = 0.947 [atm]; P0 = 1 [atm]
    
#     def Pres(P0Pa = None, AeAt = None, Pe = None, P0 = None):
        
# Purely Converging Ducts

# class Nozzle:
    
#     def mdot(state):
        
#         if isinstance(state, Isen):
#             return state.rstar * SpeedOfSound.a(state.Tstar) * state.Astar

#         raise TypeError(f"Invalid Input Type ({type(state)})")
    
    
    
#     def generate_converging(P0Pa = None, P0 = None, Pa = None, T0 = None, T = None, Ae = None):
        
#         g = config.GAMMA
#         gm = g - 1            
        
#         choked = np.exp(g/gm * np.log(gm/2 + 1))
        
#         [P0Pa, P0, Pa] = config.ratio_check(P0Pa, P0, Pa, 'Nozzle.generate_converging')
                
#         if P0Pa > choked:
#             mach1 = 1
#         else:
#             mach1 = Nozzle.Converging.from_pres(P0Pa, P0 = P0).mach
        
#         state = Nozzle.Converging(mach1, P0 = P0)
            
#         if T0 is None: state.T = T
#         else: state.T0 = T0
        
#         state.A = Ae
        
#         return state
    
#     class Converging(Isen):
        
#         def __init__(self, M, **kwargs):
#             super().__init__(M, **kwargs)
        
#         @property
#         def mdot(self):
            
#             g = config.GAMMA
#             R = config.R
            
#             return np.sqrt(g/R/self.T0) * (5/6)**3 * self.P0 * self.Astar
        
        
#         def Ax(self, x):
            
#             if (x < self.Astar):
#                 AAs_x = self.A / x
            
#             else:
#                 AAs_x = x / self.A
            
#             inter = InterIsen.from_area1(self, AAs_x, sup = False)
#             return inter     
        
    
#     class ShockPos(Enum):
#         subsonic = 1        # Pcrit1 < Pa/P0
#         shock_in_nozzle = 2 # Pcrit2 < Pa/P0 < Pcrit1
#         over_expanded = 3   # Pcrit3 < Pa/P0 < Pcrit2
#         under_expanded = 4  #          Pa/P0 < Pcrit3
    
    
#     def generate_deLaval(
#         AeAt = None, PaP0 = None, PeP0 = None, 
#         Pe = None, Pa = None, P0 = None, 
#         Ae = None, At = None):
#         # AeAt / At + Ae
#         # Pe, Pa, P0
        
#         if PaP0 is None: PaP0 = PeP0
#         elif PeP0 is None: PeP0 = PeP0
        
#         [AeAt, Ae, At] = config.ratio_check(AeAt, Ae, At, 'generate_deLaval [AeAt]')
#         [PaP0, Pa, P0] = config.ratio_check(PaP0, Pa, P0, 'generate_deLaval [PaP0]', False)
#         [PeP0, Pe, P0] = config.ratio_check(PeP0, Pe, P0, 'generate_deLaval [PeP0]', False)
        
#         crit1 = Isen.from_area(AeAt, A0=At, A=Ae, sup = False)
#         crit3 = Isen.from_area(AeAt, A0=At, A=Ae, sup = True)
#         crit1.P = Pa        
        
#         # if PaP0 is not None and config.Q_(crit1.P_P0).to('').m < PaP0:
#         #     crit1.P0_P = 1 / PaP0
        
#         inter1 = InterIsen(crit1, crit3)
#         inter2 = Normal(crit3)  
        
#         nozzle = Nozzle.deLaval(AeAt, PaP0, PeP0, inter1, inter2,
#                                 Ae = Ae, At = At,
#                                 Pa = Pa, P0 = P0, Pe = Pe)
        
#         return nozzle
    
    
    
    
#     class deLaval:
        
#         def __init__(self, AeAt, PaP0, PeP0, inter13, inter32, **kwargs):
            
#             self.AeAt = AeAt
#             self.PaP0 = PaP0
#             self.PeP0 = PeP0
            
#             self.interIsen = inter13
#             self.interNormal = inter32
            
#             self.crit1 = inter13.state1
#             self.crit3 = inter13.state2    
#             self.crit2 = inter32.state2
        
        
        
#             self.Ae = kwargs.get("Ae", None)
#             self.At = kwargs.get("At", None)
            
#             self.Pe = kwargs.get("Pe", None)
#             self.Pa = kwargs.get("Pa", None)
#             self.P0 = kwargs.get("P0", None)
            
            
#             self.area_sys = optimize(lambda M: AA0Mach(M, config.GAMMA))
        
        
        
        
#         def Ax(self, AxAt, pre = False):
            
#             state = self.crit1
            
#             # state.A_A0 -> Ae/A*
#             # self.AeAt -> Ae/At
#             # A*/At = self.AeAt / state.A_A0
#             # AxAt At/A* = Ax/A* 
            
#             AxA0 = AxAt / (self.AeAt / state.A_A0)
            
#             self.area_sys.purgeCon()

#             if False:
#                 self.area_sys.addGeq(1)
#                 mach = self.area_sys.target(3, AxA0)
            
#             else:
#                 self.area_sys.addLeq(1)
#                 mach = self.area_sys.target(0.1, AxA0)
            
#             return mach
        
        
#         @property
#         def type(self):
            
#             # if (self.Pa is None or self.Pt0 is None):
#             #     return NozzleShockPosition.subsonic
            
#             PbPt0 = self.PaP0
                    
#             if (PbPt0 > self.P1P0):
#                 return Nozzle.ShockPos.subsonic
            
#             elif (self.P2P0 < PbPt0 < self.P1P0):
#                 return Nozzle.ShockPos.shock_in_nozzle
            
#             elif (self.P3P0 < PbPt0 < self.P2P0):
#                 return Nozzle.ShockPos.over_expanded
            
#             elif (PbPt0 < self.P3P0):
#                 return Nozzle.ShockPos.under_expanded
        
#         @property
#         def Me(self):
            
#             # if (self.Pa is None or self.Pt0 is None):
#             #     return None
            
#             g = config.GAMMA
#             gm = g - 1
#             gp = g + 1
            
#             PbPt0 = self.PaP0
#             AeAt = self.AeAt
            
#             mach1 = np.sqrt(2/gm * (np.pow(PbPt0, -gm/g) - 1))
            
#             mach2_sq = -1/gm + np.sqrt(1/gm/gm + 2/gm * np.pow(2/gp, gp/gm) / np.pow(AeAt * PbPt0, 2))
#             mach2 = np.sqrt(mach2_sq)
            
#             mach3 = self.crit3.mach
#             mach4 = self.crit3.mach
            
#             shock_type = self.type
            
#             if shock_type == Nozzle.ShockPos.subsonic:
#                 return mach1
            
#             if shock_type == Nozzle.ShockPos.shock_in_nozzle:
#                 return mach2
            
#             if shock_type == Nozzle.ShockPos.over_expanded:
#                 return mach3
            
#             if shock_type == Nozzle.ShockPos.under_expanded:
#                 return mach4
        
        
#         @property
#         def mdot(self):
            
#             g = config.GAMMA
#             P0 = None
#             T0 = None
#             Astar = None
                    
#             match self.type:
#                 case Nozzle.ShockPos.subsonic:
#                     P0    = self.crit1.P0
#                     T0    = self.crit1.T0
#                     Astar = self.crit1.Astar
                
#                 case _:
#                     raise TypeError(f"Unsupported shock type ({self.type})")
            
#             if (P0 is None or T0 is None):
                
#                 if (T0 is None):
#                     T0 = self.Tt
                
#                 if (P0 is None):
#                     P0 = self.Pt0            
                
#                 if (P0 is None or T0 is None):
#                     raise ValueError(f"invalid throat value")
            
#             mdot = np.sqrt(g / config.R / T0) * (5/6)**3 * P0 * Astar
            
#             return mdot
        
#         @property
#         def P1P0(self):
#             return self.crit1.P_P0
        
#         @property
#         def P2P0(self):
#             return self.crit3.P_P0 * self.interNormal.P2_P1
        
#         @property
#         def P3P0(self):
#             return self.crit3.P_P0
"""




def isEnum(var, enumClass):
    return hasattr(var, 'name') and (var.name in enumClass.__members__)

def P02_P01_Mshock(Mshock, GAMMA = config):
    if isinstance(GAMMA, Config):
        GAMMA = GAMMA.GAMMA
    
    Ms2 = Mshock * Mshock
    return ((6*Ms2) / (Ms2 + 5))**(7/2) * (6 / (7*Ms2 - 1))**(5/2)




class Nozzle:
    
    
    def mdot(state):
        
        if isinstance(state, Isen):
            return state.rstar * SpeedOfSound.a(state.Tstar) * state.Astar
        
        raise TypeError(f"Invalid input type ({type(state)})")
        
        
    class Converging(Isen):
        
        def __init__(self, M, Ae = None, **kwargs):
            super().__init__(M, **kwargs)
            
            self.A = Ae
            self.Pa = None
            object.__setattr__(self, 'choked', False)
    
        
        @classmethod
        def shift(cls, other, Ae2, sup = False):
            
            Ae1 = other.A
            Astar = other.Astar
            
            if ( Astar > Ae2 ):
                Ae2_Ae1 = Ae1 / Ae2

                state2 = Nozzle.Converging.from_area(Ae2_Ae1, A = Ae1, sup = sup)       
                state2.A0 = Ae2
                
            else:
                A2_A0 = (Ae2/Ae1) / other.A0_A
            
                state2 = Nozzle.Converging.from_area(A2_A0, A = Ae2, sup = sup)
                state2.A0 = Astar
            
            return state2
        
        @classmethod
        def from_amb(cls, Ae, P0Pa = None, P0 = None, Pa = None, **kwargs):
    
            [P0Pa, P0, Pa] = config.ratio_check(P0Pa, P0, Pa, 'Converging.from_amb')
                
            g = config.GAMMA
            gm = g - 1
            
            if P0Pa > np.pow(1 + gm/2, g/gm): # choked
                # Aexit = A*, Pe = P*
                state = cls(1, Ae, P0 = P0)
                state.choked = True
            
            else: # not choked
                state = cls.from_pres(P0Pa, P0 = P0)
                state.A = Ae
                state.choked = False
                       
            return state
        
        def A2(self, A2):
            
            if (A2 < self.A):
                raise ValueError("Area is smaller than exit. No longer a converging Nozzle")            
            
            AA2 = A2 / self.Astar
            
            state = Nozzle.Converging.from_area(AA2, A = A2, sup = False)
            state.P0 = self.P0
            state.T0 = self.T0
                        
            return state
        
        
        def __getattr__(self, name):
            
            match name:
                
                case "Pe": 
                    return self.P
                
                case "Ae":
                    return self.A
                                
            try:
                return super().__getattr__(name)
            
            except AttributeError:
                raise AttributeError(f"Converging Nozzle does not have the attribute ({name})")    
            
        def __setattr__(self, name, value):
            
            
            match name:
                
                case "Pe":
                    self.P = value
                    return
                
                case "Pa":
                    return
        
            try:
                super().__setattr__(name, value)
            
            except AttributeError:
                raise AttributeError(f"Converging Nozzle does not have the attribute ({name})")     
        
        def __str__(self):
            
            return super().__str__() + f'choked: {self.choked}\tmdot: {(self.mdot.to_base_units()) if (self.mdot is not None) else None}\n'
    
    class ShockType(Enum):
        subsonic = 1        # Pcrit1 < Pa/P0            # Pe = Pa
        shock_in_nozzle = 2 # Pcrit2 < Pa/P0 < Pcrit1   # Pe = Pa
        over_expanded = 3   # Pcrit3 < Pa/P0 < Pcrit2   # A* = At
        under_expanded = 4  #          Pa/P0 < Pcrit3   # A* = At 
        
    class deLaval:
        
        def __init__(self, AeAt = None, Ae = None, At = None, P0Pa = None, PaP0 = None, Pa = None, P0 = None, T0 = None):
            
            if (P0Pa is None and PaP0 is not None):
                P0Pa = 1 / PaP0            
            
            [AeAt, Ae, At] = config.ratio_check(AeAt, Ae, At, 'deLaval.init')
            [P0Pa, P0, Pa] = config.ratio_check(P0Pa, P0, Pa, '', False)
            
            crit1 = Isen.from_area(AeAt, A0 = At, A = Ae, sup = False)
            crit3 = Isen.from_area(AeAt, A0 = At, A = Ae, sup = True)
            crit1.P = Pa
            
            if P0Pa is not None and config.Q_(crit1.P0_P).to('').m > P0Pa:
                crit1.P0_P = P0Pa
            
            
            self.AeAt = AeAt
            self.P0Pa = P0Pa
            
            self.interIsen = InterIsen(crit1, crit3)
            self.interNormal = Normal(crit3)
            
            self.crit1 = self.interIsen.state1
            self.crit3 = self.interIsen.state2
            self.crit2 = self.interNormal.state2
            
            self.Ae = Ae
            self.At = At
            
            self.Pa = Pa
            self.P0 = P0
            
            self.T0 = T0
            
            self.area_sys = optimize(lambda M: AA0Mach(M, config.GAMMA))
            self.mshock_sys = optimize(lambda M: P02_P01_Mshock(M, config.GAMMA))
            self.mshock_sys.addGeq(1)
            
        @property
        def type(self):
            
            if (self.P0Pa is not None):
                PaPt0 = 1 / self.P0Pa
                
                if (PaPt0 >= self.P1P0):
                    return Nozzle.ShockType.subsonic
                
                elif (self.P2P0 < PaPt0 <= self.P1P0):
                    return Nozzle.ShockType.shock_in_nozzle
                
                elif (self.P3P0 < PaPt0 <= self.P2P0):
                    return Nozzle.ShockType.over_expanded
                
                elif (PaPt0 < self.P3P0):
                    return Nozzle.ShockType.under_expanded
            
            return self._type
        
        @type.setter
        def type(self, newType):
            self._type = newType
                
                
                
            
            
        @property
        def P1P0(self): return self.crit1.P_P0
        
        @property
        def P2P0(self): return self.crit3.P_P0 * self.interNormal.P2_P1
        
        @property
        def P3P0(self): return self.crit3.P_P0
        
        
        @property
        def Me(self):
            
            g = config.GAMMA
            gm = g - 1
            gp = g + 1
            
            PaPt0 = 1 / self.P0Pa
            AeAt = self.AeAt
            
            mach1 = np.sqrt(2/gm * (np.pow(PaPt0, -gm/g) - 1))
            
            mach2_sq = -1/gm + np.sqrt(1/gm/gm + 2/gm * np.pow(2/gp, gp/gm) / np.pow(AeAt * PaPt0, 2))
            mach2 = np.sqrt(mach2_sq)
            
            mach3 = self.crit3.mach
            mach4 = self.crit3.mach
        
            
            match self.type:
                
                case Nozzle.ShockType.subsonic:
                    return mach1
                
                case Nozzle.ShockType.shock_in_nozzle:
                    return mach2
                
                case Nozzle.ShockType.over_expanded:
                    return mach3
                    
                case Nozzle.ShockType.under_expanded:
                    return mach4
                
            raise TypeError("Unknown Shock Type, cannot evaluate")
        
        
        @property
        def Mt(self):
            
            if self.AeAt < AA0Mach(self.Me):
                AtA0 = 1 / self.AeAt * AA0Mach(self.Me)
                
                return Isen.from_area(AtA0, sup = False).mach
            
            return None
              
        
        
        
        @property
        def state_exit(self):
            return Isen(self.Me, P = self.Pe)
        
        @property
        def state_throat(self):
            
            Pt0 = None
            
            if (self.Pe is not None):
                Pt0 = self.state_exit.P0
            
            return Isen(self.Mt, P0 = Pt0)
        
        
        
        @property
        def Astar(self):
            
            if not isEnum(self.type, Nozzle.ShockType):
                raise ValueError("No defined shock type")
                
            match self.type:
                
                case Nozzle.ShockType.subsonic:
                    Astar = Isen(self.Me, A = self.Ae).Astar
                
                case Nozzle.ShockType.shock_in_nozzle:
                    Astar = self.Astar2
                
                case _:
                    Astar = self.Astar1
            
            return Astar
        
        @property
        def P02_P01(self):
            
            if self.type is not Nozzle.ShockType.shock_in_nozzle:
                return 1
            
            Me2 = self.Me * self.Me
            return (1 + Me2/5)**(7/2) / self.P0Pa
        
        # AeAt * P02_P01
        @property
        def AeAs2(self):
            return self.AeAt * self.P02_P01
        
        @property
        def As2At(self):
            return 1 / self.P02_P01
        
        @property
        def Ashock_At(self):
            return Isen(self.Mshock).A_A0
        
        
        @property
        def Mshock(self):
            return self.mshock_sys.target(3, self.P02_P01)
            
        @property
        def Mshock_Normal(self):
            return Normal(Isen(self.Mshock))      
        
        
        
        
        @property
        def mdot(self):
            
            g = config.GAMMA
            T0 = self.T0
            P0 = self.P0
            Astar = self.Astar
            
            mdot = np.sqrt(g / config.R / T0) * (5/6)**3 * P0 * Astar
            
            return mdot
        

        def Ax(self, Ax, rat = False, after_throat = False, after_shock = False):
            
            if not isEnum(self.type, Nozzle.ShockType):
                raise ValueError("No defined shock type")
            
            
            if rat:
                
                match self.type:
                
                    case Nozzle.ShockType.subsonic:
                                                
                        state2 = Isen.from_area(Ax, sup = False)
                    
                    case Nozzle.ShockType.shock_in_nozzle:
                        
                        if (not after_throat) and (not after_shock): # before throat
                            state2 = Isen.from_area(Ax, sup = False)
                        
                        elif after_shock: # after throat before shock
                            state2 = Isen.from_area(Ax, sup = False)                        
                        
                        else: # after shock
                            state2 = Isen.from_area(Ax, sup = True)
                    
                    case _: # over/under expanded
                        
                        if not after_throat: # before throat
                            state2 = Isen.from_area(Ax, sup = False)
                            
                        else:
                            state2 = Isen.from_area(Ax, sup = True)
                            
                return state2.mach
                
            match self.type:
                
                case Nozzle.ShockType.subsonic:
                    
                    state = Isen.from_pres(self.P0Pa)
                    state.A = self.Ae
                    
                    Astar = state.Astar
                    
                    state2 = Isen.from_area(A = Ax, A0 = Astar, sup = False)
                
                case Nozzle.ShockType.shock_in_nozzle:
                    
                    if (not after_throat) and (not after_shock): # before throat
                        state2 = Isen.from_area(A = Ax, A0 = self.Astar1, sup = False)
                    
                    elif after_shock: # after throat before shock
                        state2 = Isen.from_area(A = Ax, A0 = self.Astar2, sup = False)                        
                    
                    else: # after shock
                        state2 = Isen.from_area(A = Ax, A0 = self.Astar1, sup = True)
                
                case _: # over/under expanded
                    
                    if not after_throat: # before throat
                        state2 = Isen.from_area(A = Ax, A0 = self.At, sup = False)
                        
                    else:
                        state2 = Isen.from_area(A = Ax, A0 = self.At, sup = True)
                        
            return state2.mach

            
            # state = self.crit1
            
            # AxA0 = AxAt / (self.AeAt / state.A_A0)
            
            # if pre:
            #     self.area_sys.addGeq(1)
            #     mach = self.area_sys.target(3, AxA0)
                
            # else:
            #     self.area_sys.addLeq(1)
            #     mach = self.area_sys.target(0.1, AxA0)
                
            # return mach
        
        
        def __getattr__(self, name):
            
            match name:
                        
                case "Pe":
                    
                    match self.type:
                        
                        case Nozzle.ShockType.over_expanded:
                            if getattr(self, 'P0', None) is None: return None
                            return self.crit3.P_P0 * self.P0
                            
                        case _:
                            if getattr(self, 'Pa', None) is None: return None
                            return self.Pa
                        
                case "Pt":
                    return self.state_throat.P    
                
                
                case "Astar1":
                    return self.At
                
                case "Astar2":
                    return self.At / self.P02_P01
                
                
                
                
                
                
                
                
                
                
                
                
                  
                        
                        
        def __str__(self):
            
            line1 = f'Msup: {self.crit1.mach}\tMsub: {self.crit3.mach}'
            line2 = f'subsonic\t[Pa/P0 > {self.P1P0}]'
            line3 = f'shock in nozzle\t[{self.P1P0} > Pa/P0 > {self.P2P0}]'
            line4 = f'over expanded\t[{self.P2P0} > Pa/P0 > {self.P3P0}]'
            line5 = f'under expanded\t[{self.P3P0} > Pa/P0]'
            
            return '\n'.join([line1, line2, line3, line4, line5])
            
            
                        
                    
        
        
        
        
        
        
        
        
        
        
        
        