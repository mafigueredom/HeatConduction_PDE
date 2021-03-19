"""
Created on Sat Jul 18 11:09:20 2020
@author: SantiagoOrtiz
"""

import numpy as np
#from ufl.operators import exp
from numpy import exp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score


def rho(T):
    """UO2 rho, Kg/m³"""
    rho_ref = 10980  # kg/m³ - reference rho at 273 K
    poly = (9.9734e-1+9.802e-6*T-2.705e-10*pow(T, 2)+4.391e-13*pow(T, 3))
    rho_T = rho_ref/pow(poly, 3)
    return rho_T


C1 = 302.27  # J/(kg*K)
C2 = 8.463e-3  # J/(kg*K²)
C3 = 8.741e7  # J/Kg
teta = 548.68  # K
Ea = 18531.7  # K


def Cp(T):
    """UO2 Cp, J/(Kg*K)"""
    part1 = pow(exp(teta/T)-1, 2)
    part2 = 2*C2*T+C3*Ea*exp(-Ea/T)/pow(T, 2)
    Cp_T = C1*pow(teta/T, 2)*exp(teta/T)/part1 + part2
    return Cp_T


def Klbda_UO2(T):
    """Thermal conductivity of
       100% dense solid UO2, W/(m k))"""
    t = T/1000
    part = np.exp(-16.35/t)*7410.5/pow(t, 5/2)
    k = 115.8/(7.5408+17.692*t+3.6142*pow(t, 2)) + part
    return k


def Klbda_ZrO2(T):
    """Thermal conductivity of
       100% dense solid ZrO2, W/(m k))"""
    k = 1.9599-T*(2.41e-4-T*(6.43e-7-T*1.946e-10))
    return k


def Klbda_UO2_poly(T, *argv):
    a, b, c, d, e, f, g = argv
    Kpoly = a*pow(T, 6)+b*pow(T, 5)+c*pow(T, 4)+d*pow(T, 3)+e*pow(T, 2)+f*T+g
    return Kpoly


Tarray = np.linspace(298, 2000, 2000)
Klbda_T = Klbda_UO2(Tarray)
guess = np.array([4E-19, -4E-15, 1E-11, -3E-08, 3E-05, -0.0273, 14.558])
popt, pcov = curve_fit(f=Klbda_UO2_poly, xdata=Tarray,
                       ydata=Klbda_T, p0=guess)

K_UO2 = lambda u : Klbda_UO2_poly(u, *popt)/1e3
K_ZrO2 = lambda u : Klbda_ZrO2(u)/1e3


if __name__ == "__main__":
    K_scifit = Klbda_UO2_poly(Tarray, *popt)
    Klbda_ZrO2 = Klbda_ZrO2(Tarray)
    plt.plot(Tarray, Klbda_UO2(Tarray), '-b')
    plt.plot(Tarray, K_scifit, '-r')
    plt.plot(Tarray, Klbda_ZrO2, '-g')
    plt.xlabel('Temperature, K')
    plt.ylabel('Kinetic constant, dimensionless')
    plt.show()

    print('Scipy R^2: {}'.format(r2_score(K_scifit, Klbda_T)))
