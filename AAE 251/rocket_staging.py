"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *

def rocket_staging(m_pay = 2000, dV = [6000, 2000], Isp = [325, 375], f_inert = [0.08, 0.1]):
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

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    rocket_staging(*[float(val) for val in sys.argv[1:]])
