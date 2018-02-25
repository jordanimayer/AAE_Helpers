"""
Students can put comments in this main docstring. Such comments will be
 important for you (and your teammates in the group project) to remember
 what your definition block does.
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *

def turbofan(rho_infty = 0.31194, p_infty = 19399, Ai = 10, V_infty = 150, Ve = 200, pe = 35000, Ae = 5, betta = 7, Ve_prime = 175, pe_prime = 30000, Ae_prime = 6):
    """
    {'username':'mayer15','assignment':'Turbofan and Ramjet Engines Intermediate Mastery','course':'fall-2016-aae-251','variables':
    {'rho_infty':'freestream density, kg/m^3',
     'p_infty':'freestream pressure, Pa',
     'Ai':'intake cross-sectional area of turbojet core, m^2',
     'V_infty':'freestream velocity of airplane, m/s',
     'Ve':'exit velocity from turbojet core, m/s',
     'pe':'exit static pressure from turbojet core, Pa',
     'Ae':'exit area of turbojet core, m^2',
     'betta':'bypass ratio, nondimensional',
     'Ve_prime':'exit velocity from duct, m/s',
     'pe_prime':'exit static pressure from duct, m/s',
     'Ae_prime':'exit area of duct, m^2',
     'm_dot_fan':'mass flow rate of air into the duct, kg/s',
     'm_dot_core':'mass flow rate of air into the turbojet core, kg/s',
     'm_dot_fuel':'mass flow rate of fuel injected into turbojet core, kg/s',
     'T':'thrust of turbofan engine, N'}}
    """

    ########
    m_dot_core = 0.5*467.91
    m_dot_fan = betta*m_dot_core
    m_dot_fuel = 0.05*m_dot_core
    T_core = (m_dot_core + m_dot_fuel)*Ve - m_dot_core*V_infty + (pe - p_infty)*Ae  #thrust from turbojet core, N
    T_fan = m_dot_fan*(Ve_prime - V_infty) + (pe_prime - p_infty)*Ae_prime  #thrust from fan, N
    T = T_core + T_fan
    print(T)
    ########

    return T

# The following allows you to run the script directly from the command line
if __name__ == '__main__':
    turbofan(*[float(val) for val in sys.argv[1:]])
