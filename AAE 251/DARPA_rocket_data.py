353# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 11:23:20 2016

@author: jordan
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *
from DARPA_rocket_mass_finder import DARPA_rocket_mass_finder
#from DARPA_rocket_mass_fraction import DARPA_dV_split_optimization
from homework21 import dV_split_optimization
from DARPA_rocket_cost import DARPA_rocket_cost
from DARPA_cheapest_rocket_mass import DARPA_cheapest_rocket_mass
from pylab import *
from DARPA_total_cost import *

def DARPA_rocket_data(m_pay=2000, dV_net=8000, Prop_list = ['APCP', 
    'LOx/Kerosene'], 
    Isp_list = [300, 353], x_guess = .50):
        
    x = x_guess
    x_temp = 0
    x_temp2 = 10000
    x_list = [] 
    m0_2 = 1000000
    m0_1 = 100000
    
    while(abs(x_temp - x) > .05 or m0_2 <= 0 or m0_1 <= 0):
        if x >= 1:
            x_list.append(x)            
            x = 0.1
        if m0_2 <=0 or m0_1 <= 0:
            x_list.append(x)
        while x in x_list:
            x += 0.05
        x_list.append(x)
        figure(0)
        clf()
        figure(1)
        clf()
        dV_1 = dV_net * x
        dV_2 = dV_net * (1-x)
        
        [Prop_2, Isp_2, m_prop_2, m_inert_2, m0_2, f_inert_2, cost_2] = DARPA_cheapest_rocket_mass(m_pay, dV_2, Prop_list[1:], Isp_list[1:])
        print('\n\nProp_2 = ', Prop_2, '\nIsp_2 = ', Isp_2, '\nm_prop_2 = ', 
              m_prop_2, '\nm_inert_2 = ', m_inert_2, '\nm0_2 = ', m0_2, '\nm0_2_lb = ', m0_2*2.20462, 
              '\nf_inert_2 = ', f_inert_2, '\ncost_2 = ', cost_2, '\n\n')
        [Prop_1, Isp_1, m_prop_1, m_inert_1, m0_1, f_inert_1, cost_1] = DARPA_cheapest_rocket_mass(m0_2, dV_1, Prop_list[:1], Isp_list[:1])
        print('\n\nProp_1 = ', Prop_1, '\nIsp_1 = ', Isp_1, '\nm_prop_1 = ', 
              m_prop_1, '\nm_inert_1 = ', m_inert_1, '\nm0_1 = ', m0_1, '\nm0_1_lb = ', m0_1*2.20462,
              '\nf_inert_1 = ', f_inert_1, '\ncost_1 = ', cost_1, '\n\n')
        
        f_inerts = [f_inert_1, f_inert_2]
        Isps = [Isp_1, Isp_2]    
        x_temp = x
        #print(f_inerts)
        #print(Isps)
        #x = DARPA_dV_split_optimization(m_pay, dV_net, Isps, f_inerts)
        x = dV_split_optimization(m_pay, dV_net, Isps, f_inerts)
        #print('\nx_temp = ', x_temp, '\nx = ', x, '\n')
        #print(x_list)
        #print(x_list)
                #print('\n\nProp_2 = ', Prop_2, '\nIsp_2 = ', Isp_2, '\nm_prop_2 = ', 
                          #m_prop_2, '\nm_inert_2 = ', m_inert_2, '\nm0_2 = ', m0_2, '\nm0_2_lb = ', m0_2*2.20462, 
                          #'\nf_inert_2 = ', f_inert_2, '\ncost_2 = ', cost_2, '\n\n')
                #print('\n\nProp_1 = ', Prop_1, '\nIsp_1 = ', Isp_1, '\nm_prop_1 = ', 
                          #m_prop_1, '\nm_inert_1 = ', m_inert_1, '\nm0_1 = ', m0_1, '\nm0_1_lb = ', m0_1*2.20462,
                          #'\nf_inert_1 = ', f_inert_1, '\ncost_1 = ', cost_1, '\n\n')
        cost = cost_1 + cost_2
        print('\nStage 1: ', Prop_1, '\nStage 2: ', Prop_2, '\nFirst-Stage Mass Fraction: ', x, '\nInitial Mass (Second Stage): ', m0_2, 
              'kg\nInitial Mass (Entire Rocket): ', m0_1, 'kg\nRocket Cost: $', group(cost))
    total_cost = DARPA_total_cost(m0_1, cost)
        
    return
    
if __name__ == '__main__':
    DARPA_rocket_data(*[float(val) for val in sys.argv[1:]])