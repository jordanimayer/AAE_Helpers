"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *

def patched_conics(mu_sun = 1.327120*10**11, AU = 1.49597870700*10**8):
    """
    {'username':'mayer15','assignment':'Patched Conics Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'mu_sun':'gravitational parameter of Sun, km^3/s^2',
     'AU':'1 Astronomical Unit, km',
     'v_infty':'hyperbolic excess speed from Earth needed to escape solar system, km/s'}}
    """

    ########
    v_earth = sqrt(mu_sun/AU)
    v_infty = (sqrt(2) - 1)*v_earth
    ########

    return v_infty

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    patched_conics(*[float(val) for val in sys.argv[1:]])
