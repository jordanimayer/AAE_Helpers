# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 21:03:16 2016

AAE 25100 - HW 14: Range and Endurance Assignment

@author: jordan
"""

from scipy import *
import sys
from scipy import integrate
from numpy import *
from pylab import *

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
    unt='ENG'  #tempoerary for when calling standard_atmo in English units   
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

def lift_coeff(W, h, V_infty, S):
    """
    Calculates lift coefficient.
    W: weight, lb
    h: altitude, ft
    V_infty: freestream velocity, ft/s
    S: reference area, ft^2
    CL: lift coefficient, dimensionless
    """
    [T, rho_infty, p] = standard_atmo(h)    #call standard_atmo to get density (rho_infty), slugs/ft^3
    q_infty = 0.5*rho_infty*V_infty**2      
    
    CL = W/(q_infty*S)
    
    return CL

def induced_drag_coeff(CL, e, S, b):
    """
    Calculates induced drag coefficient.
    CL: lift coefficient, dimensionless
    e: Oswald efficiency factor, dimensionless
    S: reference area, ft^2
    b: span, ft
    CDi: induced drag coefficient, dimensionless
    """
    AR = b**2/S     #aspect ratio, dimensionless
    
    CDi = CL**2/(pi*e*AR)
    
    return CDi
    
def turbojet_endurance(h, V_infty, TSFC, CD_0, W0, Wf, b, AR, e):
    """
    Calculates turboject endurance.
    h: altitude, ft
    V_infty: freestream velocity, ft/s
    TSFC: Thrust Specific Fuel Coefficient, lb/(lb*hr) (INCONSISTENT UNITS)
    CD_0: profile drag coefficient, dimensionless
    W0: initial weight, lb
    Wf: weight of fuel, lb
    b: span, ft
    AR: aspect ratio, dimensionless
    e: Oswald efficiency factor, dimensionless
    E: endurance, hr
    """
    S = b**2/AR     #reference area, ft^2
    CL = lift_coeff(W0, h, V_infty, S)  #lift coefficient, dimensionless
    CDi = induced_drag_coeff(CL, e, S, b)   #induced drag coefficient, dimensionless
    CD = CD_0 + CDi     #drag coefficient, dimensionless
    LoD = CL/CD     #lift-to-drag ratio, dimensionless
    Ct = TSFC / 3600    #consistent-unit TSFC, lb/(lb*s)
    W1 = W0 - Wf    #final weight, lb
    
    E = (1/Ct) * LoD * log(W0/W1)/3600
    
    return E
    
def turbojet_range(h, V_infty, TSFC, CD_0, W0, Wf, b, AR, e):
    S = b**2/AR     #reference area, ft^2
    [T, rho_infty, p] = standard_atmo(h)    #call standard_atmo to get density (rho_infty), kg/m^3
    CL = lift_coeff(W0, h, V_infty, S)  #lift coefficient, dimensionless
    CDi = induced_drag_coeff(CL, e, S, b)   #induced drag coefficient, dimensionless
    CD = CD_0 + CDi     #drag coefficient, dimensionless
    Ct = TSFC / 3600    #consistent-unit TSFC, lb/(lb*s)
    W1 = W0 - Wf    #final weight, lb
    
    #print('rho_infty*S = ', rho_infty*S, 'CL = ', CL, 'W0 = ', W0, 'W1 = ', W1)      
    
    R = 2*(2/(rho_infty*S))**0.5 * (1/Ct) * (CL**0.5/CD) * (W0**0.5 - W1**0.5)/5280
    
    return R

def turbojet_max_range(h, TSFC, CD_0, e, W0, Wf, b, AR):
    """
    Calculates turbojet range.    
    h: altitude, ft
    TSFC: Thrust Specific Fuel Coefficient, lb/(lb*hr) (INCONSISTENT UNITS)
    CD_0: profile drag coefficient, dimensionless
    W0: initial weight, lb
    Wf: weight of fuel, lb
    b: span, ft
    AR: aspect ratio, dimensionless
    e: Oswald efficiency factor, dimensionless
    R_max: maximum range, mi
    """
    S = b**2/AR     #reference area, ft^2
    Ct = TSFC / 3600    #consistent-unit TSFC, lb/(lb*s)
    W1 = W0 - Wf    #final weight, lb
    [T, rho_infty, p] = standard_atmo(h)    #call standard_atmo to find density (rho_infty), slugs/ft^3
    mLoD = ((1/3)*CD_0*pi*e*AR)**0.25/((4/3)*CD_0)        #maximum value for CL**(1/2)/CD  
    
    R_max = 2*(2/(rho_infty*S))**0.5 * (1/Ct) * mLoD*(W0**0.5-W1**0.5)/5280
    
    return R_max

def F16_range_and_endurance(h=0, W0=33500, Wf0=7000, Wf_tanks=5000, W_tanks=4800, b=32.67, AR=3, CD_0=0.015, TSFC=0.76, e=0.8):
    """
    Assesses range and endurance for F-16 Fighting Falcon (completes assignment).
    h: altitude, ft
    W0: initial weight, lb
    Wf0: weight of fuel without tanks, lb
    Wf_tanks: weight of fuel in drop tanks, lb
    W_tanks: weight of drop tanks, lb
    b: span, ft
    AR: aspect ratio, dimensionless
    CD_0: profile drag, dimensionless
    TSFC: Thrust Specific Fuel Coefficient, lb/(lb*hr) (INCONSISTENT UNITS)
    e: Oswald efficiency factor, dimensionless
    R: range, mi
    """
    #plot endurance vs. velocity at sea-level at W0
    Wf = Wf0 + Wf_tanks     #total fuel weight, lb
    V_infty_list = []   #list to store V_infty values
    E_list = []     #list to store endurances
    for V_infty in range(50, 1200):
        E = turbojet_endurance(h, V_infty, TSFC, CD_0, W0, Wf, b, AR, e)
        V_infty_list.append(V_infty)
        E_list.append(E)
    plt.figure(0)
    plot(V_infty_list, E_list)
    xlabel('Velocity, ft/s')
    ylabel('Endurance, hr')
    title('Endurance Analysis for F-16 Fighting Falcon')
    
    #Find V_infty to maximize E at sea-level at W0
    V_infty_E_max = V_infty_list[E_list.index(max(E_list))]
    
    #Find E_max at sea-level if tanks not dropped
    E_max = max(E_list)
    
    #plot range vs. velocity at sea-level at W0
    R_list = []
    for V_infty in range(50, 1200):
        R = turbojet_range(h, V_infty, TSFC, CD_0, W0, Wf, b, AR, e)
        R_list.append(R)
    plt.figure(1)
    plot(V_infty_list, R_list)
    xlabel('Velocity, ft/s')
    ylabel('Range, mi')
    title('Range Analysis for F-16 Fighting Falcon')
    
    #Find V_infty to maximize R at sea-level at W0
    V_infty_R_max = V_infty_list[R_list.index(max(R_list))]
    
    #Find R_max at sea-level if tanks not dropped
    R_max = turbojet_max_range(h, TSFC, CD_0, e, W0, Wf, b, AR)
    
    #IF fuel from tanks consumed, then tanks immediately dropped: find R_max at sea-level
    R_max_1 = turbojet_max_range(h, TSFC, CD_0, e, W0, Wf_tanks, b, AR)     #max range for period while tanks are in use
    R_max_2 = turbojet_max_range(h, TSFC, CD_0, e, W0-Wf_tanks-W_tanks, Wf0, b, AR)     #max range for period after dropping tanks
    R_max_drop = R_max_1 + R_max_2
    
    print('\nV_infty_E_max = ', V_infty_E_max, '\nE_max = ', E_max, '\nV_infty_R_max = ', V_infty_R_max, '\nR_max = ', R_max, '\nR_max_drop = ', R_max_drop)
    
    return

if __name__ == '__main__':
    F16_range_and_endurance(*[float(val) for val in sys.argv[1:]])