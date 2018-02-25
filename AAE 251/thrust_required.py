"""
This code is written by Jordan Mayer, AAE 251-002 for Homework 11. This is NOT the official code format i.e. Intermediate Masteries.
"""
from scipy import *
import sys
import pylab

def thrust_required(V_infty = linspace(50, 650, num=601), rho_infty = 2.3769*10**-3, W = 25000, CD_0 = 0.02, e = 0.8, S = 1000, b = 90):
    """
    {'username':'mayer15','assignment':'Thrust Required for Level, Unaccelerated Flight Assignment','course':'fall-2016-aae-251','variables':
    {'V_infty':'freestream velocity, ft/s',
     'rho_infty':'freestream density, slugs/ft^3',
     'W':'aircraft weight, lb',
     'CD_0':'zero-lift drag coefficient of whole aircraft, nondimensional',
     'e':'span efficiency factor, nondimensional',
     'S':'wing reference area, ft',
     'b':'span, ft',
     'TR':'thrust required, lb'}}
    """

    ########

    #Begin by calculating TR, a list of values corresponding to the list of V_infty values    

    q_infty = [(0.5*rho_infty*v**2) for v in V_infty]
    CL = [(W/(q*S)) for q in q_infty]   #lift coefficient, dimensionless
    AR = b**2/S             #aspect ratio, dimensionless
    CD_i = [(l**2/(pi*e*AR)) for l in CL]         #induced drag coefficient, dimensionless
    CD = [(CD_0+d) for d in CD_i]             #drag coefficient, dimensionless
    LoD = divide(CL, CD)        #lift/drag ratio, dimensionless    

    TR = [(W/(o)) for o in LoD]      #required thrust, lb
    
    #Plot TR vs. V_infty
    
    pylab.plot(V_infty, TR, label='Thrust Required')
    pylab.xlabel('Velocity, ft/s')
    pylab.ylabel('Thrust, lb')
    pylab.title('Aircraft Performance for NASA Aircraft During Cruise', y=1.05)
       
    #Find minimum TR and corresponding V_infty
    
    TR_min=min(TR)
    V_TR_min=V_infty[TR.index(TR_min)]
    print('TR_min = ')
    print(TR_min)
    print('V_TR_min = ')
    print(V_TR_min)
    
    #Calculate available thrust using given equation
    
    TA = [-0.01*v**2-v+5000 for v in V_infty]   #Thrust available, lb
    
    #Plot available thrust on same plot and add legend
    pylab.plot(V_infty, TA, label='Thrust Available')
    pylab.legend(loc='upper right', fancybox=True, shadow=True, bbox_to_anchor=[1.3, 1.08])
    
    #Determine range of velocities by examining when thrust required and thrust available intersect    
    
    diff = subtract(TA, TR)     #difference between thrust required and thrust available  
    V_min = 7000
    V_max = 0    
    for i in range(0, len(diff)-1, 1):
        if diff[i] >= 0.0 and diff[i-1] < 0.0:
            V_min = V_infty[i]
        elif diff[i] >= 0.0 and diff[i+1] < 0.0:
            V_max = V_infty[i]
    print('V_min = ')
    print(V_min)
    print('V_max = ')
    print(V_max)
    ########

    return TR

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    thrust_required(*[float(val) for val in sys.argv[1:]])
