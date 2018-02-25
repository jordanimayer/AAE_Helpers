# -*- coding: utf-8 -*-
"""
AAE 352 - HW 09

@author: mayer15
"""

def fatigueLife(sigma_m, sigma_a = 800, sigma_ult = 1757, sigma_f = 1937, b = -0.0762):
    """
    Calculate fatigue life with cyclic loading using Basquin's Law    
    
    sigma_m = mean stress (MPa)
    sigma_a = stress amplitude (MPa)
    sigma_ult = ultimate stress (MPa)
    sigma_f = Basquin's Law stress term (MPa)
    b = Basquin's Law exponent
    Nf = fatigue life (cycles)
    """
    
    Nf = 0.5*(sigma_a/sigma_f * 1/(1 - sigma_m/sigma_ult))**(1/b)
    print("sigma_m =", sigma_m, "\t-->\tNf =", Nf)


# Problem 2
fatigueLife(0, 800, 1757, 1937, -0.0762)
fatigueLife(200, 800, 1757, 1937, -0.0762)
fatigueLife(-200, 800, 1757, 1937, -0.0762)

# Print line to separate problems 2 and 3
print()

# Problem 3
fatigueLife(730, 210, 1233, 2030, -0.104)
fatigueLife(210, 210, 1233, 2030, -0.104)