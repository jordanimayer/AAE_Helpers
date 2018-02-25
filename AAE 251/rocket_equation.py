"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *

def rocket_equation(m_pay = 3353, dV = .753*8014, Isp = 300, f_inert = 0.08):
    """
    {'username':'mayer15','assignment':'The Rocket Equation Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'m_pay':'payload mass, kg',
     'dV':'change in velocity provided by rocket, m/s',
     'Isp':'specific impulse, s',
     'f_inert':'inert mass fraction, nondimensional',
     'm0':'initial mass of rocket, kg'}}
    """

    ########
    m0 = (m_pay*exp(dV/(9.81*Isp))*(1-f_inert))/(1-f_inert*exp(dV/(9.81*Isp)))
    ########

    print(m0)
    return m0

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    rocket_equation(*[float(val) for val in sys.argv[1:]])
