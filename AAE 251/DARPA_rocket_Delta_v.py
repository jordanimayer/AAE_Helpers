# -*- coding: utf-8 -*-
"""
AAE 251: Design Project - Rocket

This program is intended to calculate the worst-case scenario delta-V that a 
rocket must provide. It then plots this delta-V based on the altitude of the 
rocket's launch.

Created on Wed Nov 23 16:38:42 2016

@author: jordan
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *
from pylab import *

def DARPA_rocket_dV(h_launch=52976.816/3.281, h_orbit = 400, r_earth = 6378, mu_earth = 398600):
    """
    Calculates the worst-case scenario Delta-v for a rocket attempting to reach
    a circular orbit. This assumes no beneficial or harmful velocity from
    Earth (a polar orbit), and a stationary launch point. It also takes into
    account Delta-v losses from steering, drag, and gravity.
    
    h_launch: the altitude from which the rocket is launching, m
    h_orbit: the altitude of the CIRCULAR orbit in which the rocket will be, km
    r_earth: radius of the Earth, assuming the Earth is spherical in shape and
        mass distribution, km
    mu_earth: the gravitational constant of the Earth, km^3/s^2
    dV: the worst-case scenario Delta-v, km/s
    """
    
    r_orbit = h_orbit + r_earth     #the radius of the final orbit, km
    Vc = sqrt(mu_earth/r_orbit)     #the final orbital velocity, km/s
    dV_steer = .200     #dV loss from steering, km/s
    if h_launch < 20000:                
        dV_drag = 150-0.0075*h_launch       #dV loss from drag, m/s
        dV_grav = 1500-0.075*h_launch       #dV loss from gravity, m/s
    else:
        dV_drag = 0
        dV_grav = 0
    #convert dV_drag and dV_grav to km/s
    dV_drag = dV_drag/1000
    dV_grav = dV_grav/1000
    
    dV = Vc + dV_steer + dV_drag + dV_grav
    print(dV)
    
    return dV

def dV_plot():
    """
    Plots worst-case scenario Delta-v against launch altitude.
    
    h_launch: the altitude from which the rocket is launching, m
    dV: the worst-case scenario Delta-v, km/s
    """
    
    h_launch_list = linspace(40000*.304, 100000*.304)     #list of launch altitudes, ranging
                                             #from 40,000 ft to 100,000 ft
    dV_list = []
    for h in h_launch_list:
        dV_list.append(DARPA_rocket_dV(h))  #call rocket_dV to find dV 
                                            #(convert h to meters)
    plot(h_launch_list, dV_list)
    xlabel('Launch Altitude, m')
    ylabel('Delta-v Required, km/s')
    #title('Worst-Case Scenario Delta-v to Reach Orbit')
    print('dV_max = ', max(dV_list), '\ndV_min = ', min(dV_list))
    
    return
   
# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    DARPA_rocket_dV(*[float(val) for val in sys.argv[1:]])
