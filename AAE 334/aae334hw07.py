# -*- coding: utf-8 -*-
"""
AAE 334 - HW 07
Problem 4
Jordan Mayer
"""
def testSection(gamma):
    """
    M = mach number in tunnel
    gamma = specific heat ratio
    """
    To = 500 # Stagnation temperature, K
    po = 300 # Stagnation pressure, kPa
    M = 2 # Mach number
    
    # Calculate test section temperature
    T = To/(1+((gamma-1)/2)*M**2)
    
    # Calculate test section pressure
    p = po*(T/To)**(gamma/(gamma-1))
    
    print(p, " kPa\n", T, " K\n")

testSection(1.667)
testSection(1.400)
testSection(1.289)
testSection(1.044)
testSection(1.024)
