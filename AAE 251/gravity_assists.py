"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *

def gravity_assists(mu_earth = 3.986*10**5, delta = 110, v_infty = 2.0):
    """
    {'username':'mayer15','assignment':'Gravity Assist Trajectories Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'mu_earth':'gravitational parameter of Earth, km^3/s^2',
     'delta':'turn angle, deg',
     'v_infty':'hyperbolic excess speed of Earth flyby, km/s',
     'rp':'required periapsis radius to complete flyby, km'}}
    """

    ########
    a = -mu_earth/v_infty**2        #semi-major axis length, km
    print(a)
    delta = delta*pi/180        #convert delta to radians
    e = 1/(sin(delta/2))        #eccentricity
    print(e)
    rp = a*(1-e)
    print(rp)
    ########

    return rp

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    gravity_assists(*[float(val) for val in sys.argv[1:]])
