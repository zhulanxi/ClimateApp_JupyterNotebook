import numpy as np
import matplotlib.pyplot as plt
from models import taus, simpleRC, single_atten, Erad, Econv

def plot_TP(p0, T0, n, ga, a, F1, F2, Fi, kuv, kop, kir, g, taurcest):
    """
    Plots the corresponding TP profile given user input parameters.
    
    Parameters
    input (strings):
    p0: air pressure at reference level [bar]
    T0: temperature at reference level [K]
    n: pressure-IR optical depth scaling parameter, unitless
    ga: ratio of specific heats (cp/cv) 
        depending on air constituents, unitless
    a: average ratio of true lapse rate vs. 
        dry adiabatic lapse rate, unitless
    Fi: internal flux [W/m^2]
    F1(F2): TOA absorbed stellar flux in short-wave channel 1(2) [W/m^2]
    kuv: extinction coefficient in short-wave channel 1 
        (UV as in Earth's case) [m^2/kg]
    kop: extinction coefficient in short-wave channel 2 
        (optical as in Earth's case) [m^2/kg]
    kir: gray infrared extinction coefficient [m^2/kg]
    g: planetary gravitationnal acceleration, 
        used for optical depth estimation [m/s^2]
    taurcest: roughly estimated optical depth at rc boundary, unitless
    
    output:
        None    
    """
    #Check for most basic valid input
    for item in locals().items():
        try:
            float(item[1])
        except ValueError:
            print(item[0]+' must be a number.')
            return
        
    #avoiding division by 0 error
    kir = float(kir)
    if kir <= 0.:
        print(r'$\kappa_\mathrm{ir}$ must be strictly positive.')
        return
    
    p0 = float(p0)
    T0 = float(T0)
    n = float(n)
    ga = float(ga)
    a = float(a)
    F1 = float(F1)
    F2 = float(F2)
    Fi = float(Fi)
    k1 = float(kuv)/kir
    k2 = float(kop)/kir
    g = float(g)
    taurcest = float(taurcest)
    
    #Handling inputs that don't make sense
    tmp = True
    
    if p0 < 0.:
        print('p_0 must be non-negative.')
        tmp = False
    if T0 < 0.:
        print('T_0 must be non-negative.')
        tmp = False
    if (n>2. or n<1.):
        print('n takes value between 1 and 2.')
        tmp = False
    if (ga>2. or ga<1.):
        print('γ takes value between 1 and 2.')
        tmp = False
    if (a>1. or a<0.):
        print('α takes value between 0 and 1.')
        tmp = False
    if (F1<0. or F2<0. or Fi<0):
        print('All fluxes must be non-negative.')
        tmp = False
    if (k1<0. or k2<0.):
        print('All κ must be non-negative.')
        tmp = False
    if g<0.:
        print('g must be non-negative.')
        tmp = False
    if taurcest<0.:
        print('Est. τ must be non-negative.')
        tmp = False
              
    if not tmp:
        return
            
            
    
    fig, ax1 = plt.subplots(figsize=(8,6))
    ax2 = ax1.twinx()#for y-axes on both sides
        
    if k1 == 0.0000 and k2 == 0.0000:
        tau0, taurc = simpleRC(p0, T0, n, ga, a, F1+F2, Fi, kir, g, taurcest)   
    elif (k1 != 0.0000 and k2 == 0.0000):
        #Assumption: users may choose either uv or optical as the single shortwave-channel
        tau0, taurc = single_atten(p0, T0, n, ga, a, F1+F2, Fi, k1, kir, g, taurcest)
    elif (k1 == 0.0000 and k2 != 0.0000):
        tau0, taurc = single_atten(p0, T0, n, ga, a, F1+F2, Fi, k2, kir, g, taurcest)
    else:    
        tau0, taurc = taus(p0, T0, n, ga, a, F1, F2, Fi, k1, k2, kir, g, taurcest)
        
    if tau0 == -9999:
        return
    
    ps = np.linspace(0.001,p0,10000)
    #the lower limit is set to p0 by default
    #possible to become user-defined
    
    #P=P0 e^(-z/H), where H needs information about the gas to compute
    #Potential possibility to make right axis altitude
    #zs = -np.log(ps/p0)*8500
    
    ts = tau0*(ps/p0)**n
    
    prc=p0*(taurc/tau0)**(1/n)
    
    #Calculates temperature profile
    Ts=[]
    for p in ps:
        if p <= prc:
            Ts.append(Erad(p, p0, n, F1, F2, Fi, k1, k2, tau0))
        else:
            Ts.append(Econv(p, p0, T0, n, ga, a))
    
    #Calculates pressure level corresponding to the minimum temperature
    pmin = ps[Ts.index(min(Ts))]
    #tau_tropopause
    tautp = (pmin/p0)**n*tau0
    
    ax1.semilogy(Ts, ps,'r--')
    #ax2.plot(Ts, zs,linewidth=0)
    ax2.semilogy(Ts, ts, linewidth=0)

    ax1.plot([min(Ts)-20, max(Ts)+20],[prc,prc],'k:',label="RC boundary")
    ax1.plot([min(Ts)-20, max(Ts)+20],[pmin,pmin],linestyle='dashdot',color='k',label="Tropopause")

    plt.title('Atmospheric temperature-pressure profile\n Calculated parameters:'+
              r'$\tau_0=$'+
              str(round(tau0,2))+
              r', $\tau_\mathrm{rc}=$'+
              str(round(taurc,2))+'\n'+
              r', $\tau_\mathrm{tp}=$'+
              str(round(tautp,2))+'bar'+
              r', $p_\mathrm{tp}=$'+
              str(round(pmin,2)))
    
    plt.xlim(min(Ts)-20, max(Ts)+20)
    ax1.set_ylim(p0,.001)
    ax2.set_ylim(max(ts),min(ts))
    #ax2.set_ylim(min(zs),max(zs))
    ax1.set_xlabel("Temperature [K]",fontsize=12)
    ax1.set_ylabel("Pressure [bar]",fontsize=12)
    ax2.set_ylabel("Thermal optical depth",fontsize=12)
    #ax2.set_ylabel("Altitude (m)",fontsize=12)

    #Formatting
    plt.tick_params(top=True, right=True, direction='in', which='both')
    ax1.legend(fontsize=12)

    plt.show()