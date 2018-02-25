# -*- coding: utf-8 -*-
"""
AAE 251: HW 21 - Rocket Staging

Created on Mon Oct 24 14:58:04 2016

@author: mayer15
"""

from scipy import *
import sys
from scipy import integrate
from numpy import *
from pylab import *
from DARPA_rocket_mass import DARPA_rocket_mass
from DARPA_rocket_staging import DARPA_rocket_staging

def DARPA_dV_split_optimization(m_pay=2000, dV=8000, Isp=[242, 451], f_inert=[0.056, 0.136]):
    """
    This function creates a plot of the initial mass of a two-stage rocket vs.
    x, where x is the fraction of dV provided by the first stage of the rocket.
    It then finds the optimal value for x (in terms of minimizing the initial
    mass of the rocket). Two plots are created: one for 0.1 <= x <= 0.9, and a
    second with more data ponits where 0 <= x <= 1.
    'm_pay':'payload mass, kg',
    'dV':'change in velocity provided by rocket, m/s',
    'Isp':'list of specific impulses of each stage [stage 1, stage 2], s',
    'f_inert':'list inert mass fractions of each stage [stage 1, stage 2], nondimensional'
    'x_list':'list of values for x, the fraction of dV provided by the first stage of the rocket, nondimensional'
    'm0_list':'list of values for the initial mass of the rocket, kg'
    'x_opt':'the value for x that minimizes m0, dimensionless'
    """
    
    x_list = linspace(0.1, 0.9)    #for first plot
    m0_list = []
    m_prop_2_list = []
    m0_2_list = []
    m_prop_1_list = []
    m0_1_list = []
    cost_list = []
    for x in x_list:
        dV_list = [x*dV, (1-x)*dV]  #list of change in velocity provided by each stage [stage 1, stage 2], m/s
        m0_list.append(DARPA_rocket_staging(m_pay, dV_list, Isp, f_inert))
        [m_prop_2, m_inert_2, m0_2] = DARPA_rocket_mass(m_pay, dV_list[1], Isp[1], f_inert[1])
        [m_prop_1, m_inert_1, m0_1] = DARPA_rocket_mass(m_pay, dV_list[0], Isp[0], f_inert[0])
        cost = 20*(m_prop_2+m_prop_1) + 1000*m_inert_1 + 1000*m_inert_2
        m_prop_2_list.append(m_prop_2)
        m0_2_list.append(m0_2)
        m0_1_list.append(m0_1)
        cost_list.append(cost)
    
    #print(cost_list)
    #print(x_list)
    neg_ind = []
    begend = True
    for i in range(0, len(cost_list)):
        if cost_list[i] < 0:
            neg_ind.append(i)
    for i in range(0, len(neg_ind)-1):
        if neg_ind[i+1] - neg_ind[i] != 1:
            begend = False
            x_list = x_list[neg_ind[i]+1:neg_ind[i+1]]
            cost_list = cost_list[neg_ind[i]+1:neg_ind[i+1]]
            break
            """            
            x_list1 = x_list[:i]
            x_list2 = x_list[(i+1):]
            cost_list1 = cost_list[:i]
            cost_list2 = cost_list[(i+1):]
            for c in cost_list1:
                if c < 0:
                    cost_list = cost_list2
                    x_list = x_list2
                    break
            else:
                cost_list = cost_list1
                x_list = x_list1
            """
    if begend and len(neg_ind) > 0:
        if neg_ind[0] == 0:
            x_list = x_list[(neg_ind[len(neg_ind)-1]+1):]
            cost_list = cost_list[(neg_ind[len(neg_ind)-1]+1):]
        else:
            x_list = x_list[:neg_ind[0]]
            cost_list = cost_list[:neg_ind[0]]
    #print(cost_list)
    #print(x_list)
    figure(0)
    plot(x_list, cost_list)
    #xlim([0.1, 0.9])
    #ylim([0.0, 50000000])
    #ylim([0, 100000])
    xlabel('First Stage Mass Fraction (x), dimensionless')
    ylabel('Total Cost, USD')
    title('Two-Stage Rocket Cost')
    
    x_opt = x_list[cost_list.index(min(cost_list))]
    
    print('\nx_opt = ', x_opt)
    """    
    x_list2 = linspace(0, 1, 3000)
    m0_list2 = []
    for x in x_list2:
        dV_list = [x*dV, (1-x)*dV]
        m0_list2.append(rocket_staging(m_pay, dV_list, Isp, f_inert))
    figure(1)
    plot(x_list2, m0_list2)
    xlim([0, 1])
    xlabel('First-stage Fraction (x), dimensionless')
    ylabel('Initial Mass (m0), kg')
    title('Two-Stage Rocket Performance (Wider Range)')  
    """
    
    return x_opt

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    DARPA_dV_split_optimization(*[float(val) for val in sys.argv[1:]])