# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 21:38:57 2016

AAE 251-002: HW 16 - Propellers: Problem 4, parts c) and d)

@author: jordan
"""

from scipy import *
import sys
from scipy import integrate
from numpy import *
from pylab import *

def propellers(r = 1.5, omega = 120, V_infty = 100, alfa = 15):
    """
    {'username':'mayer15','assignment':'Propellers and Reciprocating Engines Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'r':'radial location of airfoil in propeller, m',
     'omega':'angular speed of propeller, rad/s',
     'V_infty':'freestream velocity of airplane, m/s',
     'alfa':'angle of attack of airfoil in propeller, deg',
     'betta':'pitch angle of airfoil in propeller, deg'}}
    """

    ########
    betta = arctan2(V_infty, (r*omega))*180/pi + alfa
    ########

    return betta

def opt_pitch_angle():
    """
    Creates a plot of radial position vs. thrust-optimal pitch angle
    for the given information.
    r_list = a list of radial position values, m
    betta_list = list of thrust-optimal pitch angle values, deg
    """
    
    r_list = linspace(0.5, 2.0, num=100)
    betta_list_1 = []   #for part c
    betta_list_2 = []   #for part d
    for r in r_list:
        betta_list_1.append(propellers(r, 120, 100, 15))
        betta_list_2.append(propellers(r, 120, 55, 15))
    plot(betta_list_1, r_list, 'k--', label='100 m/s Freestream Velocity')
    plot(betta_list_2, r_list, 'k-', label='55 m/s Freestream Velocity')
    xlabel('Thrust-Optimal Pitch Angle, deg')
    ylabel('Radial Position, m')
    title('Propeller Propulsion')
    legend(fancybox=True, shadow=True, bbox_to_anchor=[1.10, 1.00])
    
    return

if __name__ == '__main__':
    opt_pitch_angle(*[float(val) for val in sys.argv[1:]])