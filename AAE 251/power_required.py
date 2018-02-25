"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
import pylab

def power_required(V_infty, rho_infty, W, CD_0, e, S, b):
    """
    {'username':'mayer15','assignment':'Power Required for Level, Unaccelerated Flight Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'V_infty':'freestream velocity, ft/s',
     'rho_infty':'freestream density, slugs/ft^3',
     'W':'aircraft weight, lb',
     'CD_0':'zero-lift drag coefficient of whole aircraft, nondimensional',
     'e':'span efficiency factor, nondimensional',
     'S':'wing area, ft^2',
     'b':'span, ft',
     'PR':'power required, W'}}
    """

    ########
    
    #Calculate PR, PA
    
    rho_0 = 2.3769*10**-3       #density at sea-level, slugs/ft^3
    AR = b**2/S     #aspect ratio, dimensionless
    q_infty = [0.5*rho_infty*v**2 for v in V_infty]      #Bernoulli's number, kg/(m*s^2)
    CL = [W/(q*S) for q in q_infty]      #lift coefficient, dimensionless
    CD_i = [l**2/(pi*e*AR) for l in CL]  #induced drag coefficient, dimensionless
    CD = [CD_0+i for i in CD_i]      #total drag coefficient, dimensionless
    LoD = divide(CL, CD)         #lift/drag ratio, dimensionless
    TR = [W/o for o in LoD]          #required thrust, lb
    PR = multiply(TR, V_infty)    #power required, lb
    TA = [-0.01*v**2-v+5000 for v in V_infty]   #thrust available at sea level, lb
    PA = (rho_infty/rho_0) * multiply(TA, V_infty)      #power available at sea level, lb
    
    #Plot PR, PA    
    
    pylab.plt.figure(0)
    pylab.plot(V_infty, PR, label='Power Required')
    pylab.plot(V_infty, PA, '--r', label='Power Available')
    pylab.xlabel('Velocity, ft/s')
    pylab.ylabel('Power, lb')
    pylab.title('Performance Analysis for NASA Aircraft at Sea-Level')
    pylab.legend(loc='upper right', fancybox=True, shadow=True, bbox_to_anchor=[1.37, 1.00])
    
    #Find range of V_infty where flight is possible
    
    diff = subtract(PA, PR)     #difference between power required and power available  
    V_min = 7000
    V_max = 0    
    for i in range(0, len(diff)-1, 1):
        if diff[i] >= 0.0 and diff[i-1] < 0.0:
            V_min = V_infty[i]
        elif diff[i] >= 0.0 and diff[i+1] < 0.0:
            V_max = V_infty[i]
    #print('V_min = ')
    #print(V_min)
    #print('V_max = ')
    #print(V_max)
    
    ########

    return PR, PA

# The following allows you to run the script directly from the command line
#if __name__ == '__main__':
 #   power_required(*[float(val) for val in sys.argv[1:]])
