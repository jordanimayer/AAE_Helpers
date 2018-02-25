"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *

def propellers(r = 1.5, omega = 100, V_infty = 75, alfa = 15):
    """
    {'username':'mayer15','assignment':'Propellers and Reciprocating Engines Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'r':'radial location of airfoil in propeller, m',
     'omega':'angular speed of propeller, rad/s',
     'V_infty':'freestream velocity of airplane, m/s',
     'alfa':'angle of attack of airfoil in propeller, deg',
     'betta':'pitch angle of airfoil in propeller, deg'}}
    """

    ########
    betta = arctan2(V_infty, (r*omega))*180/pi + alfa
    ########

    return betta

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    propellers(*[float(val) for val in sys.argv[1:]])
