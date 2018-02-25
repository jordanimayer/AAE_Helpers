# -*- coding: utf-8 -*-
"""
AAE 339 - HW 02

Jordan Mayer, mayer15@purdue.edu
"""
import math

def cp(T):
    """
    Calculate constant-pressure specific heat capacity of the gas at a given
    temperature
    
    T: temperature, K
    cp: constant-pressure specific heat capacity, J/(kg*K)
    """
    cp = 0.959 + 1.16*10**(-4)*T + 3.65*10**(-8)*T**2
    cp = cp*1000        # convert from kJ/(kg*K) to J/(kg*K)
    
    return cp

def pe(Te):
    """
    Calculate exit pressure for a given exit temperature
    
    Te: exit temperature, K    
    pe: exit pressure, Pa
    """
    
    c = 0.959*math.log(Te/1200) + 1.16*10**(-4)*(Te-1200) # intermediate variable
    c = c + 3.65*10**(-8)/2*(Te**2 - 1200**2)
    pe = 0.180*10**6*math.exp(c)
    
    return pe
    
def h(T):
    """
    Calculate static specific enthalpy for a given temperature
    
    T: temperature, K
    cp: constant-pressure specific heat capacity, J/(kg*K)
    h: static specific enthalpy, J/kg
    """
    
    h = cp(T)*T
    
    return h

def ue(Te):
    """
    Calculate exit velocity for a given exit temperature
    
    Te: exit temperature, K
    h0: stagnation specific enthalpy, J/kg
    he: static specific enthalpy at exit, J/kg
    ue: exit velocity, m/s
    """
    
    h0 = h(1200)
    he = h(Te)
    ue = (2*(h0-he))**(1/2)
    
    return ue

R = 8.314/(30*10**(-3))     # ideal gas constant for this substance

def gamma(T):
    """
    Calculate specific heat ratio for a given temperature
    
    T: temperature, K
    gamma: specific heat ratio
    """
    
    gamma = cp(T)/(cp(T)-R)
    
    return gamma
    
def a(T):
    """
    Calculate sound speed for a given temperature
    
    T: temperature, K
    gamma: specific heat ratio
    a: sound speed, m/s
    """
    
    a = (gamma(T)*R*T)**(1/2)
    
    return a

"""
Calculate and print ue, pe, and ae for three different Te values
"""

for Te in [1100, 1000, 900]:
    print("\nConditions at Te = ", Te, " K :")
    print("ue = ", ue(Te), " m/s")
    print("pe = ", pe(Te)/(10**6), " MPa")
    print("ae = ", a(Te), " m/s")