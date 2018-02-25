"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *

def launch_losses(g = 9.81):
    """
    {'username':'mayer15','assignment':'Losses During Launch Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'g':'gravitational acceleration (assume constant), m/s^2',
     'dv_gravity_loss':'gravity loss, m/s'}}
    """

    ########
    ########
    dvg = []
    t_set = range(0,90) # provides a resolution of one second
    for t in t_set:
        thetta = (90 - t)*pi/180
        dvg.append(g*cos(thetta))
    dv_gravity_loss = sum(dvg)

    return dv_gravity_loss

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    launch_losses(*[float(val) for val in sys.argv[1:]])
