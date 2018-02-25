# -*- coding: utf-8 -*-
"""
AAE 334 - HW 11

@author: mayer15
"""

# Problem 4

def c_p(p_ratio, M1=1.7, gamma=1.4):
    """
    Determine pressure coefficient in isentropic flow.
    
    p_ratio = p/p1, dimensionless (p1 = freestream pressure)
    M1 = freestream Mach number, dimensionless
    gamma = specific heat ratio, dimensionless    
    """
    return 2/(gamma*M1**2)*(p_ratio - 1)
    
print("Cp2\t=\t", c_p(1.052))
print("Cp3\t=\t", c_p(0.7685))
print("Cp4\t=\t", c_p(1.420))
print("Cp5\t=\t", c_p(1.043))