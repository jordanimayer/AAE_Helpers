"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys

def standard_atmo(h = 2999):
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
    unt = input('What unit system are you using? (please enter "SI" if using SI units (h is in meters) or "ENG" if using English units (h is in feet)\n')
    #convert h to SI units (m) if necessary
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
    print('\n', 'T = ', T, '\n', 'p = ', p, '\n', 'rho = ', rho)
    ########

    return T, rho, p

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    standard_atmo(*[float(val) for val in sys.argv[1:]])
