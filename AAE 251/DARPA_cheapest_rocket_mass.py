# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 15:02:22 2016

@author: jordan
"""

from scipy import *
import sys
from scipy import integrate
from numpy import *
from DARPA_rocket_mass_finder import DARPA_rocket_mass_finder
from DARPA_rocket_cost import DARPA_rocket_cost

def DARPA_cheapest_rocket_mass(m_pay=2000, dV=1200, Prop_list = ['Solid', 
    'Liquid1', 'Liquid2', 'Liquid3'], Isp_list=[268, 353, 391, 339]):
    m_prop_list = []
    m_inert_list = []
    m0_list = []
    f_inert_list = []
    cost_list = []
    
    for Isp in Isp_list:
        
        Prop = Prop_list[Isp_list.index(Isp)]
        [m_prop, m_inert, m0, f_inert] = DARPA_rocket_mass_finder(m_pay, dV, Isp)        
        if Isp == 242:
            m_inert_solid = m_inert
            m_inert_liquid = 0
        else:
            m_inert_solid = 0
            m_inert_liquid = m_inert
        cost = DARPA_rocket_cost(m_prop, m_inert_solid, m_inert_liquid)
        m_prop_list.append(m_prop)
        m_inert_list.append(m_inert)
        m0_list.append(m0)
        f_inert_list.append(f_inert)
        cost_list.append(cost)
        #print(Prop, Isp, m_prop, m_inert, m0, f_inert, cost)
    
    cost = min(cost_list)
    ind = cost_list.index(cost)    
    #m0 = min(m0_list)
    #ind = m0_list.index(m0)
    m_prop = m_prop_list[ind]
    m_inert = m_inert_list[ind]
    m0 = m0_list[ind]
    f_inert = f_inert_list[ind]
    Prop = Prop_list[ind]
    Isp = Isp_list[ind]
    
    #print(Prop, Isp, m_prop, m_inert, m0, f_inert, cost)
    return [Prop, Isp, m_prop, m_inert, m0, f_inert, cost]

if __name__ == '__main__':
    DARPA_cheapest_rocket_mass(*[float(val) for val in sys.argv[1:]])

