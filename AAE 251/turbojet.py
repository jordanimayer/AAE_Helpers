"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *

def turbojet(rho_infty = 0.088909, p_infty = 5529.3, Ai = 8, V_infty = 250, Ve = 600, pe = 40000, Ae = 2):
    """
    {'username':'mayer15','assignment':'Jet Thrust and Turbojet Engines Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'rho_infty':'freestream density, kg/m^3',
     'p_infty':'freestream pressure, Pa',
     'Ai':'intake cross-sectional area, m^2',
     'V_infty':'freestream velocity of airplane, m/s',
     'Ve':'exit velocity from turbojet engine, m/s',
     'pe':'exit static pressure, Pa',
     'Ae':'exit area, m^2',
     'm_dot_air':'mass flow rate of air into the turbojet engine, kg/s',
     'm_dot_fuel':'mass flow rate of fuel injected into turbojet engine, kg/s',
     'T':'thrust of turbojet engine, N'}}
    """

    ########
    m_dot_air = rho_infty*Ai*V_infty
    m_dot_fuel = 0.05*m_dot_air
    
    T = (m_dot_air + m_dot_fuel)*Ve - m_dot_air*V_infty + (pe-p_infty)*Ae
    #print('m_dot_air = ', m_dot_air, '\nm_dot_fuel = ', m_dot_fuel, '\nT = ', T)
    ########
    
    return T

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    turbojet(*[float(val) for val in sys.argv[1:]])
