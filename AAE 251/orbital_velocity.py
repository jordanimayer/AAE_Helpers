"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *

def orbital_velocity(a = 6500, r = 6778, mu_earth = 3.986*10**5):
    """
    {'username':'mayer15','assignment':'Key People and Preliminaries Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'a':'semi-major axis, km',
     'r':'distance from center of Earth, ',
     'mu_earth':'gravitational parameter of Earth, km^3/s^2',
     'energy':'specific energy of orbit, km^2/s^2',
     'KE':'specific kinetic energy, km^2/s^2',
     'PE':'specific potential energy, km^2/s^2',
     'v':'velocity of spacecraft at r, km/s'}}
    """

    ########
    v = sqrt(mu_earth/r)
    print(v)
    ########

    return v

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    orbital_velocity(*[float(val) for val in sys.argv[1:]])
