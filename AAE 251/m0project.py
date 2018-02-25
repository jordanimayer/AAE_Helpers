# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 14:09:13 2016

@author: Jonathan
"""

from scipy import *
import sys
from scipy import integrate
from numpy import *
import matplotlib.pyplot as plt

def m0project(dV=8182, f_inert=[0.073, 0.191, 0.136], Isp=[242, 242, 353]):
    x = .01 #%dv by 1st stage
    y = .01 #%dv left from 1st stage by 2nd stage
    m_pay = 2000 #payload
    m0_list = [] #list of m0 values
    m_list = []
    m_list.append(m_pay)
    g = 9.81
    counter = 0
    x_list = []
    y_list = []
    z_list = []
    w_list = []
    while (x < .95):
        mtemp_list = []
        mtemp_list.append(m_pay)
        m0temp = []
        dV_list = []
        dV1 = x * dV
        dV_list.append(dV1)
        dV_leftover = (1-x) * dV #delta V from 2nd and 3rd stages combined
        ytemplist = []
        y = .01
        while (y < .95):
            dV2 = y * dV_leftover #dV from second stage
            dV_list.append(dV2)
            dV3 = (1-y) * dV_leftover #dV from third stage
            dV_list.append(dV3)
            while (counter < 2):
                m = mtemp_list[counter] * exp(dV_list[len(dV_list) - 1 - counter] / (g * Isp[len(Isp) - counter - 1])) * (1 - f_inert[len(f_inert) - counter - 1]) / (1-f_inert[len(f_inert) - counter - 1] * exp(dV_list[len(dV_list) - 1 - counter] / (g * Isp[len(Isp) - counter - 1])))
                #mass of current stage, starts at 3rd stage then does 2nd stage before returning to outer loop
                mtemp_list.append(m)
                counter = counter + 1            
            ytemplist.append(y)
            m0temp.append(mtemp_list[len(mtemp_list) - 1])
            y = y + .001
            mtemp_list = []
            mtemp_list.append(m_pay)
            counter = 0
        m0min_upperstages = min(m0temp)
        y_optimal = ytemplist[m0temp.index(min(m0temp))]
        m = m0min_upperstages * exp(dV1 / (g * Isp[0])) * (1 - f_inert[0]) / (1 - f_inert[0] * exp(dV1 / (g * Isp[0])))
        m0 = m + m0min_upperstages
        if (m0 > 0 and m0 < 200000 and m0min_upperstages > 0):
            m0_list.append(m0)
            m_list.append(m)
            x_list.append(x)        
            y_list.append(y_optimal * (1-x))
            z_list.append((1-y_optimal) * (1-x))
            w_list.append(y_optimal * (1-x) + (1-y_optimal) * (1-x))
        x = x + .001
    
    plt.plot(x_list, m0_list, 'b')
    plt.plot(y_list, m0_list, 'g')
    plt.plot(z_list, m0_list, 'r')
    #plt.plot(w_list, m0_list, 'g')
    plt.xlabel('% value')
    plt.ylabel('mass (kg)')
    plt.title('Initial Rocket Mass vs. % of dV by Stages')
    plt.legend(['Stage 1','Stage 2', 'Stage 3'],loc='upper right')
    plt.show
    
    m0_min = m_list[m0_list.index(min(m0_list))]
    x_opt = x_list[m0_list.index(min(m0_list))]
    y_opt = y_list[m0_list.index(min(m0_list))]
    z_opt = z_list[m0_list.index(min(m0_list))]
    print('Total m0: ',min(m0_list))
    #print('First stage m0: ', m_list[m0_list.index(min(m0_list))])
    print('First stage %dv: ', x_list[m0_list.index(min(m0_list))])
    print('Second stage %dv: ', y_list[m0_list.index(min(m0_list))])
    print('Third stage % dv: ', z_list[m0_list.index(min(m0_list))])
    
    return [x_opt, y_opt, z_opt]
    
if __name__ == '__main__':
    m0project(*[float(val) for val in sys.argv[1:]])