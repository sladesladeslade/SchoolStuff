# Standard Atmosphere Class
# Slade Brooks
# spbrooks4@gmail.com
# i stole dis from class (notes)


import sys, os
sys.path.append(os.getcwd() + r"/Spring 2025")
import numpy as np
import utils.units as uu


g0 = 9.80665                                                    # gravity (m/s^2)
R = 8314.32                                                     # J/kmol-K
W0 = 28.9644                                                    # kg/kmol
ziSTD = [0., 11., 20., 32., 47., 51., 71., 84.852]              # zi for std day (km)
ziSTD = [z*1000 for z in ziSTD]                                 # zi for std day (m)
LiSTD = [-6.5, 0., 1., 2.8, 0., -2.8, -2]                       # Li for std day (K/km)
LiSTD = [L/1000 for L in LiSTD]                                 # Li for std day (K/m)
ToSTD = 288.15                                                  # T0 for std day (K)
PoSTD = 101325.                                                 # P0 for std day (Pa)
rhooSTD = 1.225                                                 # rho0 for std day (kg/m^3)
aoSTD = 340.3                                                   # a0 for std day (m/s)
# Ti at each i for std day (K)
TiSTD = [ToSTD]
for i in range(1, len(LiSTD)+1):
    TiSTD.append(LiSTD[i - 1]*(ziSTD[i] - ziSTD[i - 1]) + TiSTD[i - 1])
# Pi at each i for std day (Pa)
PiSTD = [PoSTD]
for k in range(1, len(LiSTD)+1):
    if LiSTD[k - 1] != 0:
        PiSTD.append(PiSTD[k - 1]*(TiSTD[k - 1]/TiSTD[k])**(g0*W0/(R*LiSTD[k - 1])))
    else:
        PiSTD.append(PiSTD[k - 1]*np.exp(-g0*W0*(ziSTD[k]-ziSTD[k - 1])/(R*TiSTD[k - 1])))


class stdAtmos():
    """
    Standard atmosphere 1976 base class - valid to 86 km (282,152 ft).
   
    Methods
    -------
    T(h)
        Returns atmospheric temperature (deg F) at h (ft).
    P(h)
        Returns atmospheric pressure (psf) at h (ft).
    PR(h)
        Returns pressure ratio at h (ft).
    TR(h)
        Returns temperature ratio at h (ft).
    rho(h)
        Returns density (slugs/ft^3) at h (ft).
    dR(h)
        Returns density ratio at h (ft).
    sqrtdR(h)
        Returns square root of density ratio at h (ft).
    Aspeed(h)
        Returns speed of sound (ft/s) at h (ft).
    velA(h)
        Returns speed of sound (kt) at h (ft).
    aR(h)
        Returns speed of sound ratio at h (ft).
    qMs(h)
        Returns Q/M^2 (lb/ft^2) at h (ft).
    spW(h)
        Returns specific weight (lbm/ft^3) at h (ft).
    VRkin(h)
        Returns kinematic viscosity (ft^2/s) at h (ft).
    mu(h)
        Returns dynamic viscosity (lbm/ft*s) at h (ft).
    z(h)
        Returns geopotential altitude (m) from geometric altitude (m).
    """
   
    def __init__(self, zi:list=ziSTD, Li:list=LiSTD, Ti:list=TiSTD, Pi:list=PiSTD, T0:float=ToSTD, P0:float=PoSTD,
                rho0:float=rhooSTD, a0:float=aoSTD):
        """
        Parameters
        ----------
        zi : list
            Geo-potential altitude markers (km).
        Li : list
            Temperature lapse rate at each zi (K/m).
        Ti : list
            Temperature at each zi (K).
        Pi : list
            Pressure at each zi (Pa).
        T0 : float
            Temperature at sea level (K).
        P0 : float
            Pressure at sea level (Pa).
        rho0 : float
            Density at sea level (kg/m^3).
        a0 : float
            Speed of sound at sea level (m/s).
        """
        self.zi = zi
        self.Li = Li
        self.Ti = Ti
        self.Pi = Pi
        self.T0 = T0
        self.P0 = P0
        self.rho0 = rho0
        self.a0 = a0
       
   
    def _i(self, z:float):
        """
        Parameters
        ----------
        z : float
            Geo-potential altitude (m).
       
        Returns
        -------
        i : int
            Corresponding i value.
        """
        zi = [zs for zs in self.zi if z >= zs][-1]
        i = self.zi.index(zi)
        return i, zi
   
   
    def _start(self, h:float):
        """
        Quickly startup for doing calcs.
        """
        h = h*uu.ft2m
        zh = self.z(h)
        i, zi = self._i(zh)
        return zh, zi, i
   

    def T(self, h:float):
        """
        Returns atmospheric temperature (deg F) at h (ft).
        """
        zh, zi, i = self._start(h)
        T = self.Ti[i] + self.Li[i]*(zh - zi)
        return uu.k2degF(T)


    def P(self, h:float):
        """
        Returns atmospheric pressure (psf) at h (ft).
        """
        zh, zi, i = self._start(h)
        if self.Li[i] != 0:
            P = self.Pi[i]*(self.Ti[i]/(uu.degF2k(self.T(h))))**(g0*W0/(R*self.Li[i]))
        else:
            P = self.Pi[i]*np.exp(-(g0*W0*(zh - zi))/(R*self.Ti[i]))
        return P*uu.pa2psf
   
   
    def PR(self, h:float):
        """
        Returns pressure ratio at h (ft).
        """
        return self.P(h)/self.P(0.)
   
   
    def TR(self, h:float):
        """
        Returns temperature ratio at h (ft).
        """
        return uu.degF2k(self.T(h))/self.T0
       
       
    def rho(self, h:float):
        """
        Returns density (slugs/ft^3) at h (ft).
        """
        return self.rho0*(self.PR(h)/self.TR(h))*uu.kgm32slugft3
   
   
    def dR(self, h:float):
        """
        Returns density ratio at h (ft).
        """
        return self.rho(h)/self.rho(0.)
   
   
    def sqrtdR(self, h:float):
        """
        Returns square root of density ratio at h (ft).
        """
        return np.sqrt(self.dR(h))
   
   
    def Aspeed(self, h:float):
        """
        Returns speed of sound (ft/s) at h (ft).
        """
        return self.a0*np.sqrt(self.TR(h))*uu.ms2fts
   
   
    def velA(self, h:float):
        """
        Returns speed of sound (kt) at h (ft).
        """
        return self.Aspeed(h)*uu.fts2kt
   
   
    def aR(self, h:float):
        """
        Returns speed of sound ratio at h (ft).
        """
        return self.Aspeed(h)/self.Aspeed(0.)
   
   
    def qMs(self, h:float):
        """
        Returns Q/M^2 (lb/ft^2) at h (ft).
        """
        return 1481.354*self.PR(h)
       
   
    def spW(self, h:float):
        """
        Returns specific weight (lbm/ft^3) at h (ft).
        """
        return 32.1740484*self.rho(h)
       
   
    def VRkin(self, h:float):
        """
        Returns kinematic viscosity (ft^2/s) at h (ft).
        """
        return ((0.226968*10**(-7))*(uu.degF2r(self.T(h))**1.5))/(self.rho(h)*(uu.degF2r(self.T(h))+198.73))
    
    
    def mu(self, h:float):
        """
        Returns dynamic viscosity (lbm/ft*s) at h (ft).
        """
        return self.VRkin(h)*self.rho(h)
   
   
    @staticmethod
    def z(h:float):
        """
        Parameters
        ----------
        h : float
            Geometric height (m).
           
        Returns
        -------
        z : float
            Geo-potential altitude (m).
        """
        r = 6356.577*1000   # radius of earth (m)
        z = r*h/(r + h)     # geopot height        
        return z
    

if __name__ == "__main__":
    test = stdAtmos()
    h = 5000
    print(h, "ft")
    print("T:", test.T(h))
    print("P:", test.P(h)/144)
    print("TR:", test.TR(h))
    print("PR:", test.PR(h))
    print("rho:", test.rho(h)*10**4)
    print("dR:", test.dR(h))
    print("a:", test.Aspeed(h))
    print("a kt:", test.velA(h))
    h = 10000
    print(h, "ft")
    print("T:", test.T(h))
    print("P:", test.P(h)/144)
    print("TR:", test.TR(h))
    print("PR:", test.PR(h))
    print("rho:", test.rho(h)*10**4)
    print("dR:", test.dR(h))
    print("a:", test.Aspeed(h))
    print("a kt:", test.velA(h))
    h = 30000
    print(h, "ft")
    print("T:", test.T(h))
    print("P:", test.P(h)/144)
    print("TR:", test.TR(h))
    print("PR:", test.PR(h))
    print("rho:", test.rho(h)*10**4)
    print("dR:", test.dR(h))
    print("a:", test.Aspeed(h))
    print("a kt:", test.velA(h))
    h = 100000
    print(h, "ft")
    print("T:", test.T(h))
    print("P:", test.P(h)/144)
    print("TR:", test.TR(h))
    print("PR:", test.PR(h))
    print("rho:", test.rho(h)*10**4)
    print("dR:", test.dR(h))
    print("a:", test.Aspeed(h))
    print("a kt:", test.velA(h))
    h = 7500.
    print(h, "ft")
    print("T:", test.T(h))
    print("P:", test.P(h))
    print("TR:", test.TR(h))
    print("PR:", test.PR(h))
    print("rho:", test.rho(h))
    print("dR:", test.dR(h))
    print("a:", test.Aspeed(h))
    print("qMs:", test.qMs(h))
    print("a kt:", test.velA(h))
    print("mu:", test.mu(h))