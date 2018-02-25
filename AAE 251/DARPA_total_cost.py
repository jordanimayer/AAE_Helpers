# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 13:38:42 2016

@author: jordan
"""

from DARPA_rocket_cost import DARPA_rocket_cost

def group(number):
    s = '%d' % number
    groups = []
    while s and s[-1].isdigit():
        groups.append(s[-3:])
        s = s[:-3]
    return s + ','.join(reversed(groups))

def DARPA_total_cost(rocket_mass=62502, rocket_cost=6986887):
    system_mass = rocket_mass/0.214
    print('System Mass :', group(system_mass), 'kg' )
    plane_inert_mass = system_mass*0.55
    plane_fuel_mass = system_mass*(1-0.214-0.55)
    
    total_cost = rocket_cost+plane_inert_mass*10000+plane_fuel_mass*1000
    print('\nTotal Cost: $', group(total_cost))
    
    return(total_cost)
    
if __name__ == '__main__':
    DARPA_total_cost(*[float(val) for val in sys.argv[1:]])