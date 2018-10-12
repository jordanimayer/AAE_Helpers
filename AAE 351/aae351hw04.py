# -*- coding: utf-8 -*-
"""
AAE 351
HW 04
Problem 6

Create plot of probability mass function for waiting times
"""

import numpy
import matplotlib.pyplot as plt

# set up tau, theta list, and x list
theta_list = [4, 1, 0.5]
tau = 80
x_list = numpy.arange(0, tau, 0.1)

for theta in theta_list:
    # calculate f(x) for each x and store in list
    f_x_list = []
    for x in x_list:
        f_x = theta/tau * (1 - x/tau)**(theta-1)
        f_x_list.append(f_x)
    
    # set up title and axis labels
    ttl = 'Jordan Mayer, HW 04 Problem 6a, Theta = ' + str(theta)
    xlab = 'x'
    ylab = 'f(x; ' + str(theta) + ', ' + str(tau) + ')'
    
    # create plot
    plt.figure()
    plt.plot(x_list, f_x_list)
    plt.title(ttl)
    plt.xlabel(xlab)
    plt.ylabel(ylab)