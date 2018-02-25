"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys

def continuity_and_momentum(V = 15, omega = 10, y = 0.5):
    """
    {'username':'mayer15','assignment':'Practical Uses of Continuity and Momentum Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'V':'translational velocity of helicopter, m/s',
     'omega':'rotational speed of rotor, rad/s',
     'y':'location along rotor, m',
     'V_infty_local':'local freestream velocity, m/s'}}
    """

    ########
    #calculate V-infty_local
    
    V_infty_local = (omega * y) + V
    
    ########

    return V_infty_local

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    continuity_and_momentum(*[float(val) for val in sys.argv[1:]])
