# -*- coding: utf-8 -*-
"""
AAE 251: HW 26 - Velocity Losses

Created on Wed Nov  9 14:53:29 2016

@author: mayer15
"""

from scipy import *
import sys
from scipy import integrate
from numpy import *
from pylab import *

def circular_orbits(mu_earth = 3.986*10**5, r_earth = 6378, h_min = 100, h_max = 3000):
    """
    This function is used to calculate various characteristics of a circular orbit.
    mu_earth: MG of the earth, km^3/s^2
    r_earth: radius of Earth, assuming that Earth is spherical in shape and mass distribution, km
    h_min: minimum altitude, km
    h_max: maximum altitude, km
    r_min: minimum radius, km
    r_max: maximum radius, km
    r_list: list of values for orbital radius, km
    v_c: list of values for orbital velocity, km/s
    ke: list of values for specific kinetic energy, km^2/s^2
    pe: list of values for specific potential energy, km^2/s^2
    E: list of values for total specific energy, km^2/s^2
    """
    
    r_min = h_min + r_earth
    r_max = h_max + r_earth
    h_list = range(h_min, h_max)
    r_list = range(r_min, r_max)
    v_c = []
    ke = []
    pe = []
    E = []
    for r in r_list:
        v_c.append(sqrt(mu_earth/r))
        ke.append(mu_earth/(2*r))
        pe.append(-mu_earth/r)
        E.append(mu_earth/(2*r) - mu_earth/r)
    
    figure(0)
    plot(h_list, v_c, label='Orbital Velocity, km/s')
    ylabel('Orbital Velocity, km/s')
    xlabel('Altitude, km')
    title('Velocities of Circular Orbits')
    figure(1)
    plot(h_list, ke, label='Specific Kinetic Energy')
    plot(h_list, pe, label='Specific Potential Energy')
    plot(h_list, E, label='Specific Total Energy')
    ylabel('Energy, km^2/s^2')
    xlabel('Altitude, km')
    legend(shadow=True, fancybox=True, bbox_to_anchor = [1.50, 1.20])
    title('Specific Energies of Circular Orbits')
    
    return
        
# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    circular_orbits(*[float(val) for val in sys.argv[1:]])