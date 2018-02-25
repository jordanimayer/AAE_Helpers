"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *

def takeoff_and_landing_distance(W = 105435*9.81, T = 760000, rho_infty = 1.23, S = 371.676, CL_max = 1.5, CL = 0.1, CD = 0.005233248, mu_r = 0.02):
    """
    {'username':'mayer15','assignment':'Takeoff and Landing Distance Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'W':'weight of aircraft during takeoff, N',
     'T':'thrust of aircraft during takeoff, N',
     'rho_infty':'freestream density, kg/m^3',
     'S':'reference area, m^2',
     'CL_max':'maximum lift coefficient, nondimensional',
     'CL':'lift coefficient during takeoff, nondimensional',
     'CD':'drag coefficient during takeoff, nondimensional',
     'mu_r':'rolling friction coefficient, nondimensional',
     's_lo':'liftoff distance, m'}}
    """

    V_lo = 1.2*((2*W)/(rho_infty*S*CL_max))**0.5    #liftoff speed, m/s
    V_l = 1.3*sqrt((2*W/(rho_infty*S*CL_max)))
    print('V_l = ', V_l)
    q_infty = 0.5*rho_infty*(0.7*V_l)**2       #Bernoulli's number
    print('q_infty = ', q_infty)
    D = q_infty*S*CD    #drag force, N
    print('D = ', D)
    L = q_infty*S*CL    #lift force, N
    print('L = ', L)
    
    s_lo = (1.44*W**2)/(9.81*rho_infty*S*CL_max*(T-(D+mu_r*(W-L))))
    s_l = (1.69*W**2)/(9.81*rho_infty*S*CL_max*(D+0.4*(W-L)))

    print('s_l = ', s_l)    
    
    return s_lo

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    takeoff_and_landing_distance(*[float(val) for val in sys.argv[1:]])
