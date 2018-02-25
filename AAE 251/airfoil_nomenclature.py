"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys

def airfoil_nomenclature(LoD = 20, rho_infty = 1.23, V_infty = 30, c = 2, D = 100, M = 1000):
    """
    {'username':'mayer15','assignment':'Airfoil Nomenclature and Force Coefficients Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'LoD':'lift to drag ratio, nondimensional',
     'rho_infty':'freestream density, kg/m^3',
     'V_infty':'freestream velocity, m/s',
     'c':'chord length, m',
     'D':'drag force, N',
     'M':'aerodynamic moment, N*m',
     'cl':'lift coefficient, nondimensional',
     'cd':'drag coefficient, nondimensional',
     'cm':'moment coefficient, nondimensional'}}
    """

    ########
    q_infty = 0.5 * rho_infty * V_infty * V_infty   #Bernoulli's number, kg/(m * s^2)
    cd = D/(q_infty * c)        #drag coefficient (per unit span)
    cl = cd * LoD               #lift coefficient (per unit span)
    cm = M/(q_infty * c * c)        #moment coefficient (per unit span)
    ########

    return cl, cd, cm

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    airfoil_nomenclature(*[float(val) for val in sys.argv[1:]])
