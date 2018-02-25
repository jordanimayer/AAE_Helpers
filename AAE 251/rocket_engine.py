"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *

def rocket_engine(m_dot = 500, Ve = 3500, pe = 30000, Ae = 10):
    """
    {'username':'mayer15','assignment':'Rocket Engines and Propellants Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'m_dot':'mass flow rate of propellant, kg/s',
     'Ve':'exit velocity, m/s',
     'pe':'static pressure at exit, Pa',
     'Ae':'exit area, m^2',
     'Isp':'specific impulse, s'}}
    """

    ########
    p_infty = 1.01325*10**5      #freestream pressure at sea-level, N/m
    g_0 = 9.81      #gravitational constant at sea-level, m/s^2
    T = m_dot*Ve+(pe-p_infty)*Ae    #thrust at sea-level, N
    Isp = T/(g_0*m_dot)
    ########

    return Isp

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    rocket_engine(*[float(val) for val in sys.argv[1:]])
