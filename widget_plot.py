from ipywidgets import interactive,interact, HBox, Layout,VBox
import ipywidgets as widgets
from plotting import plot_TP

def plotFinal():
    """
    Maybe can add argument for different default planets.
    """
    interactive_plot = interactive(plot_TP, p0='1.', T0='288', n='2.', ga='1.4', a='0.6',
                                   F1='7.', F2='233.', Fi='0.0', kuv='1.6758e-2', kop='2.9792e-5', 
                                   kir='1.86e-4', g='9.8',taurcest = '0.1')
    
    #After entering value into each textbox, either click on another box
    #or hit enter to apply change.
    for i in range(13):
        interactive_plot.children[i].continuous_update=False
    
    #style1 = {'description_width': '1em'}
    style2 = {'description_width': '2em'}#2em works for most boxes
    style4 = {'description_width': '4em'}
    
    #Description does not show formatting
    for i in range(12):
        interactive_plot.children[i].style=style2
        
    interactive_plot.children[0].description=r'$p_0$: Air pressure at reference level [bar]'
    interactive_plot.children[1].description=r'$T_0$: Temperature at reference level [K]'
    interactive_plot.children[2].description=r'$n$: Pressure-IR optical depth scaling factor (unitless)'
    interactive_plot.children[3].description=r'$\gamma$: Ratio of specific heats (cp/cv) of major gas (unitless)'
    interactive_plot.children[4].description=r'$\alpha$: Average ratio of true lapse rate vs. dry adiabatic lapse rate (unitless)'
    interactive_plot.children[5].description=r'$F_1$: TOA absorbed stellar flux in short-wave channel 1 [W/m^2]'
    interactive_plot.children[6].description=r'$F_2$: TOA absorbed stellar flux in short-wave channel 2 [W/m^2]'
    interactive_plot.children[7].description=r'$F_i$: Internal flux [W/m^2]'
    interactive_plot.children[8].description=r'$\kappa_\mathrm{uv}$: Extinction coefficient in short-wave channel 1 [m^2/kg]'
    interactive_plot.children[9].description=r'$\kappa_\mathrm{op}$: Extinction coefficient in short-wave channel 2 [m^2/kg]'
    interactive_plot.children[10].description=r'$\kappa_\mathrm{ir}$: Gray infrared extinction coefficient [m^2/kg]'
    interactive_plot.children[11].description=r'$g$: Planetary gravitationnal acceleration [m/s^2]'
    
    interactive_plot.children[12].description=r'Est. $\tau_\mathrm{rc}$:Roughly estimated optical depth at RC boundary (unitless)'
    interactive_plot.children[12].style=style4

    controls = HBox(interactive_plot.children[:-1], layout = Layout(flex_flow='row wrap'))
    output = interactive_plot.children[-1]
    
    display(VBox([controls, output]))
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        