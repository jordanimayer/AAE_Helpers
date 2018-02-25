# -*- coding: utf-8 -*-
"""
AAE 251: Design Project - Rocket

Created on Fri Dec  2 21:37:57 2016

@author: jordan
"""

def DARPA_rocket_cost(m_prop, m_inert_solid, m_inert_liquid):
    """
    The purpose of this function is to calculate the total cost of the rocket 
    system.
    
    m_prop_total: the total propellant mass of all stages, kg
    m_inert_solid: the total inert mass of all solid-based stages, kg
    m_inert_liquid: the total inert mass of all liquid-based stages, kg
    rocket_cost: the total cost of the rocket system, USD
    """
    
    rocket_cost = m_prop * 20 + m_inert_solid * 500 + m_inert_liquid * 1000
    
    return rocket_cost
    # The following allows you to run the script directly from the command line
if __name__ == '__main__':
    DARPA_rocket_cost(*[float(val) for val in sys.argv[1:]])