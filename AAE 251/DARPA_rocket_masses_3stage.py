# -*- coding: utf-8 -*-
"""
AAE 251: Design Project - Rocket

Created on Fri Dec  2 21:46:51 2016

@author: jordan
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *
from pylab import *
from DARPA_cheapest_rocket_mass import DARPA_cheapest_rocket_mass
from m0project import m0project
from DARPA_rocket_staging import DARPA_rocket_staging

def DARPA_rocket_3stage(m_pay = 2000, dV_net = 8182, Prop1_list = ['APCP (Solid)', 
    'LOx/Kerosene', 'Nitrogen Tetraoxide/Hydrazine'], Prop2_list = ['NT-Hydrazine', 'LOx/Kerosene'], 
    Isp1_list = [242, 300, 292], Isp2_list = [339, 353], x_guess = .40, y_guess = .40, z_guess = .20):
    """
    The purpose of this function is to calculate masses for various stages of
    the rocket system.
    """
    x = x_guess
    y = y_guess
    z = z_guess
    x_temp = 0
    y_temp = 0
    z_temp = 0
    
    while(abs(x_temp - x) > .2*x or abs(y_temp - y) > .2*y or abs(z_temp - z) > .2*z):    
        
        if (x + y + z - 1 > 0.001):
            raise Exception('Those mass fractions do not add to 1!')
        dV_1 = dV_net * x
        dV_2 = dV_net * y
        dV_3 = dV_net * z
        
        [Prop_3, Isp_3, m_prop_3, m_inert_3, m0_3, f_inert_3, cost_3] = DARPA_cheapest_rocket_mass(m_pay, dV_3, Prop2_list, Isp2_list)
        print('\n\nProp_3 = ', Prop_3, '\nIsp_3 = ', Isp_3, '\nm_prop_3 = ', 
              m_prop_3, '\nm_inert_3 = ', m_inert_3, '\nm0_3 = ', m0_3, 
              '\nf_inert_3 = ', f_inert_3, '\ncost_3 = ', cost_3)
        [Prop_2, Isp_2, m_prop_2, m_inert_2, m0_2, f_inert_2, cost_2] = DARPA_cheapest_rocket_mass(m0_3, dV_2, Prop1_list, Isp1_list)
        print('\n\nProp_2 = ', Prop_2, '\nIsp_2 = ', Isp_2, '\nm_prop_2 = ', 
              m_prop_2, '\nm_inert_2 = ', m_inert_2, '\nm0_2 = ', m0_2, 
              '\nf_inert_2 = ', f_inert_2, '\ncost_2 = ', cost_2)
        [Prop_1, Isp_1, m_prop_1, m_inert_1, m0_1, f_inert_1, cost_1] = DARPA_cheapest_rocket_mass(m0_2, dV_1, Prop1_list, Isp1_list)
        print('\n\nProp_1 = ', Prop_1, '\nIsp_1 = ', Isp_1, '\nm_prop_1 = ', 
              m_prop_1, '\nm_inert_1 = ', m_inert_1, '\nm0_1 = ', m0_1, 
              '\nf_inert_1 = ', f_inert_1, '\ncost_1 = ', cost_1)
        
        f_inerts = [f_inert_1, f_inert_2, f_inert_3]
        Isps = [Isp_1, Isp_2, Isp_3]    
        
        x_temp = x
        y_temp = y
        z_temp = z
        [x, y, z] = m0project(dV_net, f_inerts, Isps)
        print(x, y, z)

    #print (m0_1, m0_2, m0_3)    
    
    return #[m0_3, m0_2, m0_1]

if __name__ == '__main__':
    DARPA_rocket_3stage(*[float(val) for val in sys.argv[1:]])