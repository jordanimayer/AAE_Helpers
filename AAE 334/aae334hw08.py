# -*- coding: utf-8 -*-
"""
AAE 334 - HW 08

Jordan Mayer

@author: mayer15
"""

# Problem 3

def Isp(gas, R, gamma, To = 2000, Me = 5, g = 9.81):
    """
    Find specific impulse of a rocket engine with the given characteristics.
    
    gas = name of gas (string)
    R = specific gas constant, J/(kg*K)
    gamma = specific heat ratio
    To = stagnation temperature, K
    Me = exit Mach number
    g = gravitational constant, m/s^2
    
    Isp = specific impulse, s
    """
    
    Isp = (gamma*R*To)**0.5 / g * Me / (1+(gamma-1)/2 * Me**2)**0.5
    print("Isp for", gas, ":\t", Isp, "\n")
    
Isp("hydrogen", 4124, 1.41)
Isp("air", 287, 1.4)
Isp("steam", 462, 1.33)