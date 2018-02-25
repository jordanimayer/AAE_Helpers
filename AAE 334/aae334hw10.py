# -*- coding: utf-8 -*-
"""
AAE 334 - HW 10

Jordan Mayer
@author: mayer15
"""
import math

# Problem 2

def c_p(p_ratio, M1=1.7, gamma=1.4):
    """
    Determine pressure coefficient in isentropic flow.
    
    p_ratio = p/p1, dimensionless (p1 = freestream pressure)
    M1 = freestream Mach number, dimensionless
    gamma = specific heat ratio, dimensionless    
    """
    return 2/(gamma*M1**2)*(p_ratio - 1)
    
p_rat2 = 0.9497
p_rat3 = 0.7685
p_rat4 = 1.420
p_rat5 = 1.054

# calculate pressure coefficients
c_p2 = c_p(p_rat2)
c_p3 = c_p(p_rat3)
c_p4 = c_p(p_rat4)
c_p5 = c_p(p_rat5)

print("Cp2\t=\t", c_p2)
print("Cp3\t=\t", c_p3)
print("Cp4\t=\t", c_p4)
print("Cp5\t=\t", c_p5)

# calculate normal and axial coefficients
c_n = 0.5*(c_p4 + c_p5 - c_p2 - c_p3)
eps = 3*math.pi/180    # diamond angle in radians
c_a = 0.5*(c_p2 + c_p4 - c_p3 - c_p5)*math.tan(eps)

# calculate lift and drag coefficients
alpha = 4*math.pi/180  # angle of attack in radians
c_l = c_n*math.cos(alpha) - c_a*math.sin(alpha)
c_d = c_n*math.sin(alpha) + c_a*math.cos(alpha)

# print lift and drag coefficients
print("\nc_l\t=\t", c_l, "\nc_d\t=\t", c_d)