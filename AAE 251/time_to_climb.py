"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
import pylab
from scipy.optimize import fsolve
from standard_atmo import standard_atmo
from power_required import power_required

def time_to_climb(h1 = 10000, h2 = 20000):
    """
    {'username':'mayer15','assignment':'Climbing Flight and Ceilings Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'h1':'initial altitude, ft',
     'h2':'final altitude, ft',
     'T_climb':'time to climb, min'}}
    """

    ########
    #Call power_required for PR and PA
    V_infty = linspace(50,650,601)      #velocity, ft/s
    #call standard atmosphere function
    [T0, rho0, p0] = standard_atmo(0)     #rho0 is sea-level density, slugs/ft^3
    W=25000    
    CD_0=0.02       #zero-lift drag coefficient, dimensionless
    e=0.85      #Oswald efficiency factor, demensionless
    S=1000      #wing area, ft^2
    b=95        #span, ft
    [PR,PA] = power_required(V_infty,rho0,W,CD_0,e,S,b)
   
    #Determine max RC and corresponding angle    
    
    RC = [(p/W)*60 for p in (subtract(PA,PR))]    #rate of climb, ft/min
    RC_max = max(RC)
    print('\nRC_max = ', RC_max)
    theta = arcsin(RC_max/(V_infty[RC.index(RC_max)]*60))*180/pi    #climb angle, degrees
    print('\ntheta = ', theta)
    
    #Determine service ceiling
    
    # The following provides the rate of climb for each altitude of interest in one foot intervals
    def rate_of_climb(h):
        [T, rho, p] = standard_atmo(h)  #Temperature, density, and pressure at h
        PA_alt = [a*(rho/rho0) for a in PA]    #Available power, adjusted for altitude
        PR_alt = [r*(rho0/rho)**0.5 for r in PR]   #Required power, adjusted for altitude
        R_C = [(d/W)*60 for d in (subtract(PA_alt,PR_alt))]     #set of rates of climb at h for varying speeds
        
        R_C_max = max(R_C)      #maximum rate of cilmb at H
        return R_C_max-100
    def y(x):
        y = x**2 - 2
        return y
    SC = fsolve(rate_of_climb, 30000)   #service ceiling, feet
    print('\nSC = ', SC[0])
    
    #Plot max RC vs. altitude from sea-level to service ceiling
    
    RC_mx = []      #list of max RCs
    h_set = []      #list of altitudes
    times=0
    for h in range(0, int(SC[0]), 100):
        R_C_temp = rate_of_climb(h)+100  #rate of climb, ft/min
        RC_mx.append(R_C_temp)
        h_set.append(h)
    pylab.plt.figure(1)
    pylab.plot(h_set, RC_mx)
    pylab.title('Performance Analysis for DARPA Aircraft')
    pylab.xlabel('Altitude, ft')
    pylab.ylabel('Rate of Climb, ft/min')
    
    #Determine minimum time to climb from sea-level to service ceiling
    
    One_over_RC=[1/rc for rc in RC_mx]  #one over rate of climb (min/ft)
    T_climb=integrate.simps(One_over_RC,h_set)  #time to climb, min
    print('\nFrom sea-level to service ceiling: T_climb = ', T_climb)
    
    # note, to use Simpson's integration rule, you would use this:
    # integral_value = integrate.simps(y, x)
    # where y is a list of y-values and x is a list of x-values

    #Determine range flown during glide from SC to 10,000ft

    delh = int(SC[0])-10000        #change in altitude, ft
    LoD = 12.5      #lift/drag ratio, dimensionless
    R = LoD * delh/5280  #glide range, miles
    print('\nRange: ',R)

    #Determine time to climb from 10,000 ft to service ceiling

    RC_mx = []
    h_set = []
    for h in range(10000, int(SC[0]), 100):
        R_C_temp = rate_of_climb(h)+100  #rate of climb, ft/min
        RC_mx.append(R_C_temp)
        h_set.append(h)
    One_over_RC=[1/r for r in RC_mx]
    T_climb = integrate.simps(One_over_RC,h_set)
    print('\nFrom 10000 feet to service ceiling: T_climb = ', T_climb)

    ########

    return T_climb

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    time_to_climb(*[float(val) for val in sys.argv[1:]])
