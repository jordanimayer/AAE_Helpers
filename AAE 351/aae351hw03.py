# -*- coding: utf-8 -*-
"""
Jordan Mayer
AAE 35103
HW 03

Perform basic arithmetic.
"""

import numpy as np
from scipy.special import comb as nCr
import matplotlib.pyplot as plt

"""
Question 2
"""
print('\nQuestion 2:\n')
# Part a: check PMF properties
x_list = range(1,10)
p_x_list = np.log10((x_list + np.ones(9))/x_list)
print('Max p(x) = ', max(p_x_list))
print('Min p(x) = ', min(p_x_list))
print('Sum(p(x)) = ', sum(p_x_list))

# Part b: calculate each p(x)
for x, p_x in zip(x_list, p_x_list):
    print('p(', x, ') = ', round(p_x, 4))

"""
Question 4
"""
print('\n\nQuestion 4:\n')
# Calculate P(X = x) for each x from 0 through 5
x_list = range(0,6)
p_list = []
for x in x_list:
    p = nCr(5, x) * 0.75**x * 0.25**(5-x)
    p_list.append(p)
    print(p)
plt.plot(x_list, p_list)
plt.title('Jordan Mayer: HW 03 Problem 4')
plt.xlabel('X')
plt.ylabel('PMF')
plt.show()