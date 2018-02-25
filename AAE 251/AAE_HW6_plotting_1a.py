# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 14:53:42 2016

@author: mayer15
"""

"""
AAE 251
HW 6
Problem 1
Part a

Goal: Plot C_N, C_A, C_L, C_D, and LoD, as function of angle of attack

rn: nose radius, ft
rc: base cone radius, ft
delc: cone half angle, degrees
C_N: normal coefficient, dimensionless
C_A: axial coefficient, dimensionless
C_L: lift coefficient, dimensionless
C_D: drag coefficient, dimensionless
LoD: lift/drag ratio, dimensionless
"""

import pylab
import scipy
import numpy
import math

rn = 0.3
rc = 1.25
delc_deg = 70
    
#calculate values
    
delc = delc_deg * (math.pi / 180)   #convert delc to radians
angle = numpy.linspace(0, delc, 50)   #list of angles of attack, radians
C_N = [((1 - (rn / rc)**2 * (math.cos(delc))**2) * (math.cos(delc))**2 * 
    math.sin(2 * a)) for a in angle]
C_A = [((1 - (math.sin(delc))**4) * (rn / rc)**2 + (2 * (math.sin(delc))**2 
    * (math.cos(a))**2 + (math.cos(delc))**2 * (math.sin(a))**2) * 
    (1 - (rn / rc)**2 * (math.cos(delc))**2)) for a in angle]
C_L = numpy.multiply(C_N, [math.cos(a) for a in angle]) - numpy.multiply(C_A, [math.sin(a) for a in angle])
C_D = numpy.multiply(C_N, [math.sin(a) for a in angle]) + numpy.multiply(C_A, [math.cos(a) for a in angle])
LoD = numpy.divide(C_L, C_D)
angle_deg = numpy.multiply(angle, (180/math.pi))    #angle of attack in degrees

#compose plot

pylab.plot(angle_deg, C_N, 'r-', label='Normal Coefficient')  #plot C_N as a function of alpha
pylab.plot(angle_deg, C_A, 'b-', label='Axial Coefficient')  #plot C_A as a function of alpha
pylab.plot(angle_deg, C_L, 'c-', label='Lift Coefficient')  #plot C_L as a function of alpha
pylab.plot(angle_deg, C_D, 'g-', label='Drag Coefficient')  #plot C_D as a function of alpha
pylab.plot(angle_deg, LoD, 'k-', label='Lift/Drag Ratio')  #plot LoD as a function of alpha
pylab.axis([0, delc_deg, -2.0, 2.0])
pylab.xlabel('Angle of Attack (degrees)')
pylab.ylabel('Calculated Values (dimensionless)')
pylab.title('Hypersonic Aerodynamics of Sphere-Cone Entry System')
pylab.legend(loc='upper center', bbox_to_anchor=(1.24, 1.05), fancybox=True, shadow=True)