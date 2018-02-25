"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *

def range_and_endurance(rho_infty = 1.23, S = 40, b = 20, TSFC = 0.86, CD0 = 0.02, e = 0.85, W0 = 50000, Wf = 15000):
    """
    {'username':'your_username','assignment':'Range and Endurance Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'rho_infty':'freestream density, kg/m^3',
     'S':'reference area, m^2',
     'b':'span, m',
     'TSFC':'thrust specific fuel consumption - NONCONSISTENT, N/(N*hr)',
     'CD0':'zero-lift drag coefficient of whole aircraft, nondimensional',
     'e':'span efficiency factor, nondimensional',
     'W0':'initial weight of the aircraft, N',
     'Wf':'FUEL weight initially onboard the aircraft, N',
     'R':'range of aircraft, m'}}
    """

    ########
    # Student's code and comments go here.
    ########

    return R

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    range_and_endurance(*[float(val) for val in sys.argv[1:]])
