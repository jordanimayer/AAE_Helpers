"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *

def open_orbits(e = 2, v_infty = 2.0, mu_earth = 3.986*10**5):
    """
    {'username':'mayer15','assignment':'Ground Tracks and Launch Sites Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'e':'eccentricity, nondimensional',
     'v_infty':'hyperbolic excess speed, km/s',
     'mu_earth':'gravitational parameter of Earth, km^3/s^2',
     'a':'semi-major axis, km',
     'rp':'periapsis radius, km',
     'vp':'periapsis velocity, km/s'}}
    """

    ########
    a = -mu_earth/(v_infty**2)
    rp = a*(1-e)
    vp = (v_infty**2 + 2*mu_earth/rp)**0.5
    print(a, '\n', rp, '\n', vp)
    ########

    return a, rp, vp

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    open_orbits(*[float(val) for val in sys.argv[1:]])
