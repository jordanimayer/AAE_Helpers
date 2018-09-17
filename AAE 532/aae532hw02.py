"""
Jordan Mayer
AAE 532
HW 02

Perform various mathematical operations to compute motion of astronomical
bodies.
"""

"""
Problem 1
"""

# define distances, km (planets/moon based on semi-major axes)
r_sun_earth = 149649952
r_earth_moon = 384400
r_earth_sc = 1.5*10**6
r_sun_jupiter = 779067093
r_sun_moon = r_sun_earth + r_earth_moon
r_sun_sc = r_sun_earth + r_earth_sc

G = 6.674*10**(-20)  # universal gravitational constant, N * km^2 / kg^2

# define masses, kg (based on gravitational parameters)
m_sun = 132712200000 / G
m_moon = 4902.801076 / G
m_earth = 398600.4418 / G
m_sc = 830
m_jupiter = 126686535 / G

"""
Part a
"""
# calculate distance of center of mass from sun, km
r_sun_cm = ((m_earth * r_sun_earth + m_moon * r_sun_moon + m_sc * r_sun_sc +
             m_jupiter * r_sun_jupiter) /
            (m_sun + m_earth + m_moon + m_sc + m_jupiter))
r_sun = 0 - r_sun_cm
r_earth = r_sun_earth - r_sun_cm
r_moon = r_sun_moon - r_sun_cm
r_sc = r_sun_sc - r_sun_cm
r_jupiter = r_sun_jupiter - r_sun_cm
print('r_sun_cm =', "{:,}".format(r_sun_cm), 'km\n')
print('r_sun =', "{:,}".format(r_sun), 'km')
print('r_earth =', "{:,}".format(r_earth), 'km')
print('r_moon =', "{:,}".format(r_moon), 'km')
print('r_sc =', "{:,}".format(r_sc), 'km')
print('r_jupiter =', "{:,}".format(r_jupiter), 'km')

print('\n')

"""
Part b
"""
"""
inertial_accel:
    calculate acceleration of one body due to another using inertial reference
    inputs:
        m_j: mass of "other" body (kg)
        r_j: position magnitude of "other" body (km)
        r_i: position magnitude of accelerating body (km)
"""
def inertial_accel(m_j, r_j, r_i):
    r_j_i = r_i - r_j
    return -G * m_j / r_j_i**2

# calculate acceleration due to various bodies
r_ddot_sun_sc = inertial_accel(m_sun, r_sun, r_sc)
r_ddot_earth_sc = inertial_accel(m_earth, r_earth, r_sc)
r_ddot_moon_sc = inertial_accel(m_moon, r_moon, r_sc)
r_ddot_jupiter_sc = inertial_accel(m_jupiter, r_jupiter, r_sc)
print('r_ddot_sun_sc =', r_ddot_sun_sc)
print('r_ddot_earth_sc =', r_ddot_earth_sc)
print('r_ddot_moon_sc =', r_ddot_moon_sc)
print('r_ddot_jupiter_sc =', r_ddot_jupiter_sc)

# calculate net acceleration
print('\nr_ddot_sc =', sum([r_ddot_sun_sc, r_ddot_earth_sc, r_ddot_moon_sc,
                            -r_ddot_jupiter_sc]), 'km/s^2')
    # negatives account for direction
    
print('\n')
    
"""
Problem 2
"""

"""
Part b
"""
"""
dominant_term:
    compute dominant term of acceleration of one body relative to another
    inputs:
        m_i: mass of accelerating body (kg)
        m_q: mass of "relative to" body (kg)
        r_q_i: distance between accelerating body and "relative to" body (km)
"""
def dominant_term(m_i, m_q, r_q_i):
    return -G * (m_i + m_q)/r_q_i**2

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

# compute MAGNITUDE of each term
d_earth_sc = dominant_term(m_sc, m_earth, r_sc - r_earth)
dp_sc_sun = perturbing_term(m_sun, r_sun - r_sc)
idp_earth_sun = perturbing_term(m_sun, r_sun - r_earth)
dp_sc_moon = perturbing_term(m_moon, r_moon - r_sc)
idp_earth_moon = perturbing_term(m_moon, r_moon - r_earth)
dp_sc_jupiter = perturbing_term(m_jupiter, r_jupiter - r_sc)
idp_earth_jupiter = perturbing_term(m_jupiter, r_jupiter - r_earth)
print('d_earth_sc =', d_earth_sc, 'km/s^2')
print('dp_sc_sun =', dp_sc_sun, 'km/s^2')
print('idp_earth_sun =', idp_earth_sun, 'km/s^2')
print('dp_sc_moon =', dp_sc_moon, 'km/s^2')
print('idp_earth_moon =', idp_earth_moon, 'km/s^2')
print('dp_sc_jupiter =', dp_sc_jupiter, 'km/s^2')
print('idp_earth_jupiter =', idp_earth_jupiter, 'km/s^2\n')

# print net perturbing terms (+/- to account for directions)
print('p_sun =', -dp_sc_sun + idp_earth_sun, 'km/s^2')
print('p_moon =', -dp_sc_moon - idp_earth_moon, 'km/s^2')
print('p_jupiter =', dp_sc_jupiter - idp_earth_jupiter, 'km/s^2\n')

# print net relative acceleration (+/- to account for directions)
print('r_ddot_earth_sc =', sum([d_earth_sc, -dp_sc_sun, idp_earth_sun,
                                -dp_sc_moon, -idp_earth_moon, dp_sc_jupiter,
                                -idp_earth_jupiter]), 'km/s^2')
    # signs account for direction