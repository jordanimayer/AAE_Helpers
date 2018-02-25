# -*- coding: utf-8 -*-
"""
AAE 251 HW 20: Rocket Equation

Created on Fri Oct 21 15:05:09 2016

@author: mayer15
"""

from scipy import *
import sys
from scipy import integrate
from numpy import *
from pylab import *

def SLS(m_pay = 10000):
    """
    Plots change in velocity versus initial mass of SLS and calculates the maximum change in velocity.
    m_pay: mass of payload, kg
    m0_list: list of initial mass values, kg
    Isp: specific impulse, s
    f_inert: inert mass fraction, dimensionless
    g0: gravitational constant at sea-level on Earth
    dV_list: list of change in velocity over launch values, m/s
    dV_max: maximum change in velocity over launch, m/s
    """

    m0_list = linspace(0, 3000000, 100)
    g0 = 9.81
    Isp = 375
    f_inert = 0.09

    dV_list = []
    for m in m0_list:
        dV_list.append(g0*Isp*log(m/(m_pay*(1-f_inert)+m*f_inert)))
    
    plot(m0_list, dV_list)
    title('SLS Performance: 10,000 kg Payload')
    xlabel('Initial Mass, kg')
    ylabel('Change in Velocity over Launch, m/s')
    print('dV = ', max(dV_list))
    
    return
    
if __name__ == '__main__':
    SLS(*[float(val) for val in sys.argv[1:]])