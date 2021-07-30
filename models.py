import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.special import gamma, gammaincc
from scipy.optimize import fsolve

D = 1.66
sigma = 5.67e-8

def taus(p0, T0, n, ga, a, 
         F1, F2, Fi, k1, k2, kir, 
         g, taurcest=1):
    """
    Calculates optical depth at reference level and at the rc boundary
    for two short-wave channels with attenuation.
    
    Parameters
    input (floats):
    p0: air pressure at reference level [bar]
    T0: temperature at reference level [K]
    n: pressure-IR optical depth scaling parameter, unitless
    ga: ratio of specific heats (cp/cv) 
        depending on air constituents, unitless
    a: average ratio of true lapse rate vs. 
        dry adiabatic lapse rate, unitless
    Fi: internal flux [W/m^2]
    F1(F2): TOA absorbed stellar flux in short-wave channel 1(2) [W/m^2]
    k1(k2): ratio of short-wave attenuation in channel 1(2)
            vs. thermal attenuation, unitless
    kir: gray infrared extinction coefficient [m^2/kg]
    g: planetary gravitationnal acceleration, 
        used for optical depth estimation [m/s^2]
    taurcest: roughly estimated optical depth at rc boundary, unitless
        
    output:
        optical depths at reference level & at rc boundary
    """
    #4b/n factor
    bn = 4*(a*(ga-1)/ga)/n
    
    def equations(t):
        tau0, taurc = t
        return ((F1/2.)*(1+D/k1+(1-D/k1)*sp.exp(-k1*taurc))+\
               (F2/2.)*(1+D/k2+(1-D/k2)*sp.exp(-k2*taurc))+\
               (Fi/2)*(2+D*taurc)-\
               sigma*T0**4*sp.exp(D*taurc)*(sp.exp(-D*tau0)+(1/(D*tau0)**bn)*\
               (gammaincc(1+bn,D*taurc)*gamma(1+bn)-gammaincc(1+bn,D*tau0)*gamma(1+bn))),
               (F1/2.)*(1+D/k1+(k1/D-D/k1)*sp.exp(-k1*taurc))+\
               (F2/2.)*(1+D/k2+(k2/D-D/k2)*sp.exp(-k2*taurc))+\
               (Fi/2)*(2+D*taurc)-\
               sigma*T0**4*(taurc/tau0)**bn)
    
    
    #kir: thermal extinction coefficient at p0, 
    #tau0 = kir(at p0) * p0/g
    
    tau0est = kir*p0*100000/g #100000 for bar to Pa conversion

    solutions = fsolve(equations,(tau0est,taurcest))

        
    print("Two short-wave channels with attenuation.")
    
    return (solutions)

def simpleRC(p0, T0, n, ga, a, F, Fi, kir, g, taurcest=1):
    """
    Calculates optical depth at reference level and at the rc boundary
    for the single channel without attenuation.
    
    Parameters (different from above)
    input (float):
    F: TOA absorbed net stellar flux (sum of F1 and F2) [W/m^2]
    """
    bn = 4*(a*(ga-1)/ga)/n
    
    def equations(t):
        tau0, taurc = t
        return ((F+Fi)*(1+D*taurc)/2. - sigma*T0**4*(taurc/tau0)**bn,
               sigma*T0**4*sp.exp(D*taurc)*(sp.exp(-D*tau0)+(1/(D*tau0)**bn)*\
               (gammaincc(1+bn,D*taurc)*gamma(1+bn)-gammaincc(1+bn,D*tau0)*gamma(1+bn)))-\
               (2+D*taurc)*(F+Fi)/2.)
    
    tau0est = kir*p0*100000/g
    
    solutions = fsolve(equations,(tau0est,taurcest))
    
    print("Simplest case without short-wave attenuation.")
    
    return (solutions)

def single_atten(p0, T0, n, ga, a, F, Fi, k, kir, g, taurcest=1):
    """
    Calculates optical depth at reference level and at the rc boundary
    for the single channel model with attenuation.
    
    Parameters (different from above)
    input (float):
    F: TOA absorbed net stellar flux [W/m^2]
    k: ratio of short-wave attenuation vs. thermal attenuation, unitless
    """
    bn = 4*(a*(ga-1)/ga)/n
    
    def equations(t):
        tau0, taurc = t
        return ((F/2.)*(1+D/k+(1-D/k)*sp.exp(-k*taurc))+\
               (Fi/2)*(2+D*taurc)-\
               sigma*T0**4*sp.exp(D*taurc)*(sp.exp(-D*tau0)+(1/(D*tau0)**bn)*\
               (gammaincc(1+bn,D*taurc)*gamma(1+bn)-gammaincc(1+bn,D*tau0)*gamma(1+bn))),
               (F/2.)*(1+D/k+(k/D-D/k)*sp.exp(-k*taurc))+\
               (Fi/2)*(2+D*taurc)-\
               sigma*T0**4*(taurc/tau0)**bn)
    
    tau0est = kir*p0*100000/g
    
    solutions = fsolve(equations,(tau0est,taurcest))
    
    print("Single short-wave channel with attenuation.")
    
    return (solutions)

#Calculates the T-P profile in two regions

def Erad(p, p0, n, F1, F2, Fi, k1, k2, tau0):
    
    tau = tau0*(p/p0)**n
    F = F1+F2
    
    if k1 == 0.0 and k2 == 0.0:
        return (((1+D*tau)*(F+Fi)/2.)/sigma)**(1/4)
    
    elif (k1 != 0.0 and k2 == 0.0):#single channel
        return ((((F/2)*(1+D/k1+(k1/D-D/k1)*np.exp(-k1*tau))+
                (Fi/2)*(1+D*tau))/sigma)**(1/4))
    elif (k1 == 0.0 and k2 != 0.0):#single channel
        return ((((F/2)*(1+D/k2+(k2/D-D/k2)*np.exp(-k2*tau))+
                (Fi/2)*(1+D*tau))/sigma)**(1/4))
    else:
        return ((((F1/2)*(1+D/k1+(k1/D-D/k1)*np.exp(-k1*tau))+
                (F2/2)*(1+D/k2+(k2/D-D/k2)*np.exp(-k2*tau))+
                (Fi/2)*(1+D*tau))/sigma)**(1/4))

def Econv(p, p0, T0, n, ga, a):
    return T0*(p/p0)**((a*(ga-1))/ga)









