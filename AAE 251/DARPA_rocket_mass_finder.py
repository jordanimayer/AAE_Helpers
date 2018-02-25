# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 10:32:19 2016

@author: jordan
"""
from scipy import *
import sys
from scipy import integrate
from numpy import *
from pylab import *
from DARPA_rocket_mass import DARPA_rocket_mass
from DARPA_closest_analogy import DARPA_closest_analogy

def DARPA_rocket_mass_finder(m_pay=3368, dV=.753*8014, Isp=451):
    m_pay_lbs = m_pay/0.453592  #payload mass, lbs  
    
    Isp_list = [284, 293, 204, 302, 286, 302, 291, 434, 319, 316, 337, 327, 303, 
                337,431, 327, 265, 289, 315, 451, 359, 319, 304, 292, 409, 294, 
                284,321, 321, 423, 427, 320, 310, 444, 466, 446, 278, 319, 265, 
                326,275, 331, 293, 297, 326, 259, 255, 231, 273, 247, 227, 281, 
                239, 282, 249, 230]
    m_pay_list = [71851, 1800, 10946, 27237, 48890, 241390, 49190, 68745, 90731,
                  529092, 76516, 159122, 4800, 19644, 68124, 48330, 8667, 8895,
                  8650, 20000, 29542, 31405, 1418800, 6479, 1100000, 310208,
                  68785, 44093, 44093, 89800, 35200, 940, 3968, 10200, 44093,
                  44093, 213, 6165, 1115, 6415, 98592, 11023, 330, 3090, 5974,
                  12842.9, 1521, 528, 3478.2, 2178, 200, 734.2, 2150, 100,
                  6305, 25]         #payload mass data, lbs
    m0_list = [328097, 264364, 74252, 251657, 220140, 466041, 241390, 469986, 
               303796, 1521193, 748554, 529092, 49190, 334591, 446216, 
               159122, 39757, 26167, 71581, 69794, 90731, 126079, 6257000,
               11366, 11731893, 1320160, 28285980, 68785, 68124, 353200,
               1417900, 10946, 48890, 41539, 106609, 81542, 1015, 21493, 8697,
               57276, 709272, 24912, 3428, 14639, 19644, 36325.1, 11794, 10758, 
               12842.9, 12924, 1562, 3478.2, 27572, 734.2, 61207, 681]            #initial mass data, lbs
    f_inert_list = [.038, .053, .056, .056, .056, .060, .060, .066, .068, .068,
                    .068, .070, .071, .072, .073, .073, .076, .078, .080, .085,
                    .085, .088, .093, .094, .096, .099, .104, .107, .110, .111,
                    .112, .112, .115, .118, .118, .123, .135, .136, .140, .142,
                    .144, .151, .160, .160, .177, .191, .205, .215, .218, .226,
                    .225, .241, .246, .281, .322, .451]       #inert mass fraction data
    #find initial guess for m0
    m0_lb = m0_list[m_pay_list.index(min(m_pay_list, key=lambda x:abs(x-m_pay_lbs)))]    #initial guess for m0, lb 
    m0_closest = DARPA_closest_analogy(Isp_list, m0_list, m0_lb, Isp)   #closest value to m0 in m0_list
    temp = 0   
    counter = 0
    m0_calc_list = []
    f_inert_calc_list = []
    #print(Isp)
    while (abs(m0_lb - temp) > .10*m0_lb):        
        if (counter > 20 and abs(m0_lb - temp) < .15*m0_lb):
            break
        elif (counter > 20 and abs(m0_lb - temp) < .20*m0_lb):
            break
        elif (counter > 20 and abs(m0_lb - temp) < .25*m0_lb):
            break
        elif (counter > 100):
            raise Exception('Infinite loop in DARPA_rocket_mass_finder!')
        #find corresponding f_inert
        f_inert = f_inert_list[m0_list.index(m0_closest)]   #corresponding inert mass fraction from f_inert_list
        #calculate m0 using f_inert and m_pay
        temp = m0_lb   #closest m0 with current guess
        [m_prop, m_inert, m0] = DARPA_rocket_mass(m_pay, dV, Isp, f_inert)    #new guess for m0
        m0_lb = m0 * 2.20462
        m0_closest = DARPA_closest_analogy(Isp_list, m0_list, m0_lb, Isp)   #closest m0 with new guess
        #print('\ntemp = ', temp/2.20462, '\nf_inert = ', f_inert, '\nm0 = ', m0_lb/2.20462, '\nm0_closest = ', m0_closest/2.20462, '\n\n')
                
        if (counter > 20):
            if (m0 not in m0_calc_list):
                m0_calc_list.append(m0)
                f_inert_calc_list.append(f_inert)
            else:
                f_inert = median(sort(f_inert_calc_list))
                [m_prop, m_inert, m0] = DARPA_rocket_mass(m_pay, dV, Isp, f_inert)
                break
        
        counter += 1
        #repeat process until closest m0 is found
    
    #print(m_prop, m_inert, m0, f_inert)
    return [m_prop, m_inert, m0, f_inert]
    
if __name__ == '__main__':
    DARPA_rocket_mass_finder(*[float(val) for val in sys.argv[1:]])