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

def rocket_staging(m_pay, dV, Isp, f_inert):
    """
    {'username':'mayer15','assignment':'Rocket Staging Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'m_pay':'payload mass, kg',
     'dV':'list of change in velocity provided by each stage [stage 1, stage 2], m/s',
     'Isp':'list of specific impulses of each stage [stage 1, stage 2], s',
     'f_inert':'list inert mass fractions of each stage [stage 1, stage 2], nondimensional',
     'm0':'initial mass of two-stage rocket, kg'}}
    """

    ########
    g0 = 9.81       #gravitational constant at sea-level on Earth, m/s^2
    m0_2 = (m_pay*exp(dV[1]/(g0*Isp[1]))*(1-f_inert[1]))/(1-f_inert[1]*exp(dV[1]/(g0*Isp[1])))  #initial mass of second stage
    m0 = (m0_2*exp(dV[0]/(g0*Isp[0]))*(1-f_inert[0]))/(1-f_inert[0]*exp(dV[0]/(g0*Isp[0])))     #initial mass of rocket
    ########

    return m0

def dV_split_optimization(m_pay=2000, dV=7767, Isp=[300, 353], f_inert=[0.085, 0.136]):
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
    for x in x_list:
        dV_list = [x*dV, (1-x)*dV]  #list of change in velocity provided by each stage [stage 1, stage 2], m/s
        m0_list.append(rocket_staging(m_pay, dV_list, Isp, f_inert))
    neg_ind = []
    begend = True
    for i in range(0, len(m0_list)):
        if m0_list[i] < 0:
            neg_ind.append(i)
    for i in range(0, len(neg_ind)-1):
        if neg_ind[i+1] - neg_ind[i] != 1:
            begend = False
            x_list = x_list[neg_ind[i]+1:neg_ind[i+1]]
            m0_list = m0_list[neg_ind[i]+1:neg_ind[i+1]]
            break
    if begend and len(neg_ind) > 0:
        if neg_ind[0] == 0:
            x_list = x_list[(neg_ind[len(neg_ind)-1]+1):]
            m0_list = m0_list[(neg_ind[len(neg_ind)-1]+1):]
        else:
            x_list = x_list[:neg_ind[0]]
            m0_list = m0_list[:neg_ind[0]]
    
    figure(0)
    plot(x_list, m0_list)
    xlim([0.1, 0.9])
    xlabel('First-stage Fraction (x), dimensionless')
    ylabel('Initial Mass (m0), kg')
    title('Two-Stage Rocket Performance')
    
    x_opt = x_list[m0_list.index(min(m0_list))]
    print(x_opt)
    """
    print('\nx_opt = ', x_opt, '\nm0_min = ', min(m0_list1))
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
    dV_split_optimization(*[float(val) for val in sys.argv[1:]])