"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *

def orbital_properties(a = 10000, e = 0.2, mu_earth = 3.986*10**5):
    """
    {'username':'mayer15','assignment':'Conic Sections Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'a':'semi-major axis, km',
     'e':'orbital eccentricity, nondimensional',
     'rp':'periapsis radius, km',
     'vp':'periapsis velocity, km/s',
     'ra':'apoapsis radius, km',
     'va':'apoapsis velocity, km/s'}}
    """

    ########
    rp = a*(1-e)
    ra = a*(1+e)
    vp = sqrt(-mu_earth/a + 2*mu_earth/rp)
    va = sqrt(-mu_earth/a + 2*mu_earth/ra)
    print('\nrp = ', rp, '\nra = ', ra, '\nvp = ', vp, '\nva = ', va)
    ########

    return rp, vp, ra, va

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    orbital_properties(*[float(val) for val in sys.argv[1:]])
