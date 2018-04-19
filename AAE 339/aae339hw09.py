# -*- coding: utf-8 -*-
"""
AAE 339 - HW 09
Problem 3 (rocket staging)

Created on Wed Apr 11 15:42:36 2018

@author: mayer15
"""

import math

def glow(deltaV_fractions, Isps, mass_fractions, m_payload = 5000, g = 9.81, 
         deltaV_total = 9000):
    """
    Simple function to calculate GLOW (first-stage inert mass) for a staged
    rocket.
    
    Parameters:
        deltaV_fractions: list of deltaV fractions of each stage from bottom
            to top, e.g. [0.30, 0.70]
        Isps: Isps of each stage, in same order as deltaV_fractions, m/s
        mass_fractions: mass fractions of each stage, in same order as previous
            lists
        m_payload: payload mass, kg
        g: graviational constant, m/s^2
        deltaV_total: total deltaV for entire rocket, m/s
        
    Return:
        GLOW of entire rocket, kg
    """
    # check for valid inputs
    if sum(deltaV_fractions) != 1:
        raise ValueError('deltaV fractions don\'t add up to 1!')
    if any(lam >= 1 for lam in mass_fractions):
        raise ValueError('mass fraction is greater than or equal to 1!')
    
    m_totals = 0.0;  # will hold total mass values for each stage    
    
    for dV_frac, Isp, lam in zip(reversed(deltaV_fractions), reversed(Isps),
                                 reversed(mass_fractions)):
        # reverse to go from top to bottom
        dV = dV_frac * deltaV_total
        x = math.exp(dV/(g*Isp))
        f = lam/(1-lam)
        m_in = m_payload*(x-1)/(1+f-x)
        m_p = f*m_in

        # to check
        m_o = m_in + m_p + m_payload
        m_f = m_in + m_payload
        
        m_payload = m_in + m_p + m_payload  # total mass of this stage is 
                                            # payload mass of next stage
        m_totals = [m_payload, m_totals]
    
    return m_totals[0]  # GLOW = total mass of lowest stage + upper stages

splits = [[0.30, 0.70], [0.50, 0.50], [0.70, 0.30]]
Isps = [280, 440]
lams = [0.92, 0.90]
for split in splits:
    print('\ndeltaV split (1 2):', split[0], split[1])
    gl = glow(split, Isps, lams)
    print('GLOW = %.0f kg' % gl)
