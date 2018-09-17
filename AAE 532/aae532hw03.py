# -*- coding: utf-8 -*-
"""
Jordan Mayer
AAE 532
HW 03

Perform various mathematical operations to compute motion of astronomical
bodies.
"""

import math

# define distances, km
r_sun_jupiter = 779067093
r_jupiter_asteroid = r_sun_jupiter
r_sun_asteroid = math.sqrt(r_sun_jupiter**2 + r_jupiter_asteroid**2)
    # by Pythagorean Theorem

G = 6.674*10**(-20)  # universal gravitational constant, N * km^2 / kg^2

# define masses, kg (based on gravitational parameters)
m_sun = 132712200000 / G
m_jupiter = 126686535 / G
m_asteroid = 100 / G

"""
dominant_term:
    compute dominant term of acceleration of one body relative to another
    inputs:
        m_i: mass of accelerating body (kg)
        m_q: mass of "relative to" body (kg)
        r_q_i: distance between accelerating body and "relative to" body (km)
"""
def dominant_term(m_i, m_q, r_q_i):
    return G * (m_i + m_q)/r_q_i**2

"""
perturbing_term:
    compute perturbing term of acceleration of one body relative to
    another due to perturbing body
    inputs:
        m_j: mass of perturbing body (kg)
        r_i_j: distance from body of interest to perturbing body (km)
               direct: body of interest = accelerating body
               indirect: body of interest = "relative to" body
"""
def perturbing_term(m_j, r_i_j):
    return G * m_j/r_i_j**2

# --- PROBLEM 1 ---
print('Problem 1:\n')
# Part a: compute acceleration of asteroid relative to sun
print('Part a:\n')
d = dominant_term(m_asteroid, m_sun, r_sun_asteroid)
dp = perturbing_term(m_jupiter, r_jupiter_asteroid)
  # Note: r_asteroid_jupiter = r_jupiter_asteroid (only magnitude)
ip = perturbing_term(m_jupiter, r_sun_jupiter)
print('Dominant term mag = ', d, 'km/s^2')
print('Dominant term in each direction = ', d/math.sqrt(2), 'km/s^2')
  # components should have equal magnitude by definitions of r_hat_1, r__hat_2
print('Direct perturbing term = ', dp, 'km/s^2')
print('Indirect perturbing term = ', ip, 'km/s^2')
net_pert_1 = dp  # net perturbing acc in r_hat_1 dir
net_pert_2 = ip  # net perturbing acc in r_hat_2 dir
print('Net perturbing mag = ', math.sqrt(dp**2 + ip**2), 'km/s^2')
net_acc_1 = d/math.sqrt(2) + dp  # net acc in r_hat_1 dir
net_acc_2 = d/math.sqrt(2) + ip  # net acc in r_hat_2 dir
print('Net acc in each direction = ', net_acc_1, 'km/s^2')
print('Net acc mag = ', math.sqrt(net_acc_1**2 + net_acc_2**2),
      'km/s^2')
  # since dp terms are equal and dp = ip

# Part b: repeat part a, but now relative to Jupiter
print('\nPart b:\n')
d = dominant_term(m_asteroid, m_jupiter, r_jupiter_asteroid)
dp = perturbing_term(m_sun, r_sun_asteroid)
  # Note: r_asteroid_sun = r_sun_asteroid (only magnitude)
ip = perturbing_term(m_sun, r_sun_jupiter)
  # Note: r_jupiter_sun = r_sun_jupiter (only magnitude)
print('Dominant term = ', d, 'km/s^2')
print('Direct perturbing term mag = ', dp, 'km/s^2')
print('Direct perturbing term in each direction = ', dp/math.sqrt(2), 'km/s^2')
print('Indirect perturbing term = ', ip, 'km/s^2')
net_pert_1 = -dp/math.sqrt(2)  # net perturbing acc in r_hat_1 dir
net_pert_2 = ip - dp/math.sqrt(2)  # net perturbing acc in r_hat_2 dir
print('Net perturbing 1 = ', net_pert_1, 'km/s^2')
print('Net perturbing 2 = ', net_pert_2, 'km/s^2')
print('Net perturbing mag = ', math.sqrt(net_pert_1**2 + net_pert_2**2),
      'km/s^2')
net_acc_1 = net_pert_1  # net acc in r_hat_1 dir
net_acc_2 = net_pert_2 - d  # net acc in r_hat_2 dir
print('Net acc 1 = ', net_acc_1, 'km/s^2')
print('Net acc 2 = ', net_acc_2, 'km/s^2')
print('Net acc mag = ', math.sqrt(net_acc_1**2 + net_acc_2**2), 'km/s^2')
theta_p_b = math.atan(-net_pert_1/net_pert_2) * 180/math.pi
print('Net pert acc angle = ', theta_p_b, 'deg')