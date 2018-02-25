# -*- coding: utf-8 -*-
"""
AAE 251: Design Project - Rocket

Created on Fri Dec  2 21:22:55 2016

@author: jordan
"""

from scipy import *
import sys
from scipy import integrate
from numpy import *

def DARPA_rocket_mass(m_pay = 2000 , dV = 8014, Isp = 383, f_inert = .116):
    """
    The purpose of this function is to calculate the propellant mass, inert 
    mass, and total initial mass of a rocket stage.
    m_pay: the payload mass, kg
    
    dV: the Delta-v the stage must perform, m/s
    Isp: the specific impulse of the rocket, s
    f_inert: the estimated inert mass fraction of the rocket
    g0: the gravitational constant at Earth's surface, m/s^2
    m_prop: the propellant mass, kg
    m_inert: the inert mass, kg
    m0: the initial mass, kg
    """

    ########
    g0 = 9.81
    x = dV/(g0*Isp)
    m_prop = (m_pay*(exp(x) - 1)*(1-f_inert))/(1-f_inert*exp(x))
    m_inert = (m_pay*f_inert*(exp(x) - 1))/(1 - f_inert*exp(x))
    m0 = m_pay+m_prop+m_inert
    m02 = (m_pay*exp(x)*(1-f_inert)/(1-f_inert*exp(x)))
    #print('m_prop = ', m_prop, '\nm_inert = ', m_inert, '\nm0 = ', m0, '\nm02 = ', m02)
    #print('propellant mass frac' ,(m_prop/m0))
    ########

    #print(m0)
    return [m_prop, m_inert, m0]

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    DARPA_rocket_mass(*[float(val) for val in sys.argv[1:]])