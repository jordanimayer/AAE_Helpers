# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 20:29:07 2016

@author: jordan
"""

from scipy import *
from scipy import integrate
from scipy.optimize import fsolve
import pylab
import sys

def standard_atmo(h = 0):
    """
    {'username':'mayer15','assignment':'The Standard Atmosphere Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'h':'geopotential altitude, m',
     'T':'temperature at h, K',
     'rho':'density at h, kg/m^3',
     'p':'pressure at h, Pa'}}
    """

    # Atmosphere data - do not change!
    T_set = [288.16,216.66,216.66,282.66,282.66,165.66,165.66,225.66] # list (array) of temperature points that define endpoints of each layer (starting at the ground), K
    h_set = [0,11000,25000,47000,53000,79000,90000,105000] # list (array) of altitude points that define endpoints of each layer (starting at the ground), m
    a_set = [-6.5*10**-3,0,3*10**-3,0,-4.5*10**-3,0,4*10**-3] # list (array) of gradient layer slopes (starting at the ground), K/m

    ########

    #prompt user for unit system (SI vs. Engish)
    #unt = input('What unit system are you using? (please enter "SI" if using SI units (h is in meters) or "ENG" if using English units (h is in feet)\n')
    #convert h to SI units (m) if necessary
    unt='ENG'    
    if unt == 'ENG':
        h = h / 3.28
    
    #important constants
    T0 = T_set[0]           #temperature at sea level, K
    p0 = 101325.0        #pressure at sea level, Pa
    rho0 = 1.22500             #density at sea level, kg/m^3
    g0 = 9.8                #gravitational constant, m/s^2
    R = 287.058                 #gas constant for dry air, J/(K*kg)

    #base values arrays
    p_set = []      #list of pressure points that define endpoints of each layer (starting at ground), Pa
    rho_set = []    #list of density points that define endpoints of each layer (starting at ground), kg/m^3
    for t in range(0, 8):
        if t == 0:
            p_set.append(p0)
            rho_set.append(rho0)
        elif t % 2 == 0:
            p_set.append(p_set[t-1] * exp((-g0 / (R * T_set[t])) * (h_set[t] - h_set[t-1])))
            rho_set.append(rho_set[t-1] * (p_set[t] / p_set[t-1]))
        else:
            p_set.append(p_set[t-1] * (T_set[t] / T_set[t-1])**(-g0 / (a_set[t-1] * R)))
            rho_set.append((rho_set[t-1] * p_set[t] * T_set[t-1]) / (p_set[t-1] * T_set[t]))
    
    #calculate T, rho, and p for h in metric units (m)  
    for x in h_set:
        indx = h_set.index(x)
        if indx != len(h_set) and h >= x and h < h_set[indx+1]:
            T1 = T_set[indx]
            p1 = p_set[indx]
            rho1 = rho_set[indx]
            h1 = h_set[indx]
            a = a_set[indx]
            if a == 0:
                T = T1
                p = p1 * exp((-g0 / (R * T)) * (h - h1))
                rho = rho1 * (p / p1)
            else:
                T = T1 + (h - h1) * a       #temperature, K
                p = p1 * (T / T1)**(-g0 / (a * R))     #pressure, Pa
                rho = (rho1 * p * T1) / (p1 * T)      #density, kg/m^3
    
    #convert T, rho, p to English units if necessary
    if unt == 'ENG':
        T = T * (9/5)               #convert from K to degrees R
        rho = 0.00194 * rho      #convert from kg/m^3 to slugs/ft^3
        p = 0.020885 * p           #convert from Pa to lb/ft^2
    #print('\n', T, '\n', p, '\n', rho)
    ########

    return T, rho, p

def power_required(v, rho_infty, W, CD_0, e, S, b):
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
    
    [T_0, rho_0, p_0] = standard_atmo(0)
    AR = b**2/S     #aspect ratio, dimensionless
    q_infty = 0.5*rho_infty*v**2     #Bernoulli's number, kg/(m*s^2)
    CL = W/(q_infty*S)      #lift coefficient, dimensionless
    CD_i = CL**2/(pi*e*AR)  #induced drag coefficient, dimensionless
    CD = CD_0+CD_i      #total drag coefficient, dimensionless
    LoD = CL/CD         #lift/drag ratio, dimensionless
    TR = W/LoD          #required thrust, lb
    PR = TR*v    #power required, lb
    TA = -0.01*v**2-v+5000   #thrust available at sea level, lb
    PA = (rho_infty/rho_0) * TA*v      #power available at sea level, lb
    
    #Plot PR, PA    
    
    #pylab.plot(V_infty, PR, label='Power Required')
    #pylab.plot(V_infty, PA, '--r', label='Power Available')
    #pylab.xlabel('Velocity, ft/s')
    #pylab.ylabel('Power, lb')
    #pylab.title('Performance Analysis for NASA Aircraft at Sea-Level')
    #pylab.legend(loc='upper right', fancybox=True, shadow=True, bbox_to_anchor=[1.37, 1.00])
    #pylab.show()    
    
    #Find range of V_infty where flight is possible
    
    #diff = subtract(PA, PR)     #difference between power required and power available  
    #V_min = 7000
    #V_max = 0    
    #for i in range(0, len(diff)-1, 1):
     #   if diff[i] >= 0.0 and diff[i-1] < 0.0:
      #      V_min = V_infty[i]
       # elif diff[i] >= 0.0 and diff[i+1] < 0.0:
        #    V_max = V_infty[i]
    #print('V_min = ')
    #print(V_min)
    #print('V_max = ')
    #print(V_max)
    
    ########

    return PR, PA

def rate_of_climb(h, W, S, b, e, CD_0, v):
    """
    h: altitude, ft
    W: weight, lb
    S: reference area, ft^2
    b: span, ft
    e: Oswald efficiency factor, dimensionless
    CD_0: zero-lift drag coefficient, dimensionless
    V_infty: range of speeds, ft/s
    RC_max: maximum rate of climb, ft/min
    """
    [T0, rho0, p0] = standard_atmo(0)    #atmospheric values at sea-level
    [T, rho, p] = standard_atmo(h)      #atmospheric values at altitude h
    [PR, PA] = power_required(v, rho0, W, CD_0, e, S, b)    #required and available power at sea-level 
    PA_alt = PA * (rho/rho0)   #available power at altitude h
    PR_alt = PR * (rho0/rho)**0.5  #required power at altitude h
    diff = PA_alt-PR_alt
    R_C = (diff/W) * 60  #rates of climb, ft/min
    
    return R_C

def max_RC_DARPA(h):
    """
    h: altitude, ft
    DARPA_maxRC: maximum rate of climb for DARPA aircraft at altitude h, ft/min
    """
    R_C_list = []
    for v in linspace(50,650,100):
        #print('max_RC_DARPA: ', v)        
        R_C = rate_of_climb(h, 25000, 1000, 95, 0.85, 0.02, v)
        R_C_list.append(R_C)
    DARPA_maxRC = max(R_C_list)
    #DARPA_maxRC = 0
    #for r in R_C:
    #    if r > DARPA_maxRC:
    #        DARPA_maxRC = r
    
    #for r in range(0, len(R_C)):
    #    if R_C[r] > DARPA_maxRC:
    #        DARPA_maxRC = R_C[r]
    
    return DARPA_maxRC-100
    

def time_to_climb(h1, h2):
    """
    {'username':'mayer15','assignment':'Climbing Flight and Ceilings Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'h1':'initial altitude, ft',
     'h2':'final altitude, ft',
     'T_climb':'time to climb for DARPA aircraft, min'}}
    """
    One_over_RC=[]
    Height=[]
    for h in range(h1, h2+1):
        R_C = max_RC_DARPA(h)+100  #rate of climb at altitude h, ft/min
        One_over_RC.append(1/R_C)   #one over rate of climb, min/ft
        Height.append(h)    #altitude, ft
    T_climb = integrate.simps(One_over_RC, Height)
    
    return T_climb

def glide_range(h1, h2, LoD_glide):
    """
    h1: initial altitude, ft
    h2: final altitude, ft
    LoD_glide: lift-to-drag ratio for glide, dimensionless
    R: range, mi
    """
    R = (h1-h2)*LoD_glide/5280
    
    return R

def DARPA_aircraft(W=25000, S=1000, b=95, e=0.85, CD_0=0.02, V_infty=linspace(50,650,100), LoD_glide=12.5):
    """
    W: weight, lb
    S: reference area, ft^2
    b: span, ft
    e: Oswald efficiency factor, dimensionless
    CD_0: zero-lift drag coefficient, dimensionless
    V_infty: range of speeds, ft/s
    LoD_glide: lift-to-drag ratio for gliding
    R_C_max_0: maximum rate of cilmb at sea-level, ft/min
    theta_0: climb angle at sea-level for maximum rate of climb, degrees
    SC: service ceiling, ft
    T_climb_0_SC_min: minimum time to climb from sea-level to service ceiling, min
    R_SC_10k: range flown by gliding from service ceiling to 10,000 ft, ft
    T_climb_10k_SC_min: minimum time to climb from 10,000 ft to service ceiling, min
    """
    #obtain sea-level values
    [T_0, rho_0, p_0] = standard_atmo(0)   #standard atmosphere values at sea-level
    PR_0 = []
    PA_0 = []
    for v in V_infty:
        #print('DARPA_aircraft: ', v)
        [PR_0_temp, PA_0_temp] = power_required(v, rho_0, W, CD_0, e, S, b)    #required and available power at sea-level 
        PR_0.append(PR_0_temp)
        PA_0.append(PA_0_temp)
    
    #plot PR vs V_infty and PA vs V_infty at sea-level
    pylab.plot(V_infty, PR_0, label='Power Required')
    pylab.plot(V_infty, PA_0, '--r', label='Power Available')
    pylab.xlabel('Velocity, ft/s')
    pylab.ylabel('Power, lb')
    pylab.title('Performance Analysis for DARPA Aircraft at Sea-Level')
    pylab.legend(loc='upper right', fancybox=True, shadow=True, bbox_to_anchor=[1.37, 1.00])
    pylab.show()
    
    #calculate maximum rate of climb at sea-level
    R_C = []
    for v in V_infty:   
        R_C_temp = rate_of_climb(0, W, S, b, e, CD_0, v)   #list of rates of climb for different V_infty values at sea-level
        R_C.append(R_C_temp)
    
    R_C_max_0 = max(R_C)
    
    #calculate corresponding climb angle
    theta_0 = arcsin(R_C_max_0/(V_infty[R_C.index(R_C_max_0)]*60))*180/pi
    
    #calculate service ceiling
    SC = int(fsolve(max_RC_DARPA, 30000)[0])
    
    #plot RC_max vs. h from sea-level to service ceiling
    h_list = linspace(0,SC,100)
    R_C_max_list = [max_RC_DARPA(h)+100 for h in h_list]
    pylab.plot(h_list, R_C_max_list)
    pylab.xlabel('Altitude, ft')
    pylab.ylabel('Maximum Rate of Climb, ft/min')
    pylab.title('Aircraft Performance for DARPA Aircraft')
    
    #calculate minimum time to cilmb from sea-level to service ceiling
    T_climb_0_SC_min = time_to_climb(0, SC)
    
    #calculate range flown by gliding from service ceiling to 10,000 ft
    R_SC_10k = glide_range(SC, 10000, LoD_glide)
    
    #calculate minimum time to cilmb from 10,000 ft to service ceiling
    T_climb_10k_SC_min = time_to_climb(10000, SC)
    
    print('\nR_C_max_0 = ', R_C_max_0, '\ntheta_0 = ', theta_0, '\nSC = ', SC, '\nT_climb_0_SC_min = ', T_climb_0_SC_min, '\nR_SC_10k = ', R_SC_10k, '\nT_climb_10k_SC_min = ', T_climb_10k_SC_min)    
    
    return R_C_max_0, theta_0, SC, T_climb_0_SC_min, R_SC_10k, T_climb_10k_SC_min
 
if __name__ == '__main__':
    DARPA_aircraft(*[float(val) for val in sys.argv[1:]])      
    