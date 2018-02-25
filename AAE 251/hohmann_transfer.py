"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *
from pylab import *

def hohmann_transfer(h1 = 400, h2_min=-6378, h2_max=0, r_earth = 6378, mu_earth = 3.986*10**5):
    """
    {'username':'mayer15','assignment':'Orbital Maneuvers Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'h1':'altitude of circular orbit 1, km',
     'h2':'altitude of circular orbit 2, km',
     'r_earth':'radius of circular Earth, km',
     'mu_earth':'gravitational parameter of Earth, km^3/s^2',
     'dV':'delta V of Hohmann transfer, km/s'}}
     NOTE: Function is modified for Design Project to only compute dV to enter transfer orbit (deorbit burn)
    """

    ########
    #Maneuver 1: enter transfer orbit
    r1 = h1 + r_earth
    h2_list = range(h2_min, h2_max+1)
    r2_list = [h2+r_earth for h2 in h2_list]
    dV_1_list = []
    for r2 in r2_list:
        Vc_1 = (mu_earth/r1)**0.5
        a_t = (r1 + r2)/2
        E_t = -(mu_earth/(2*a_t))
        Vt_p = (2*(E_t + mu_earth/r1))**0.5
        dV_1 = Vt_p - Vc_1
        dV_1_list.append(abs(dV_1))
    dV_1_min = min(dV_1_list)
    r2_opt = r2_list[dV_1_list.index(dV_1_min)]
    print('dV_1_min = ', dV_1_min)
    print('r2_opt = ', r2_opt)
    plot(r2_list, dV_1_list)
    xlabel('Periapsis Radius of Transfer Orbit, km')
    ylabel('dV Required to Enter Transfer Orbit, km/s')
    #title('Possible Hohmann Transfers to Deliver Payload')
    """
    #Maneuver 2: exit transfer orbit
    Vc_2 = (mu_earth/r2)**0.5
    Vt_a = (2*(E_t + mu_earth/r2))**0.5
    dV_2 = Vc_2 - Vt_a
    
    #Total delta V
    dV = dV_1 + dV_2
    """
    ########

    return dV_1

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    hohmann_transfer(*[float(val) for val in sys.argv[1:]])
