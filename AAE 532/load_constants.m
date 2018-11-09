%%%%%
% Jordan Mayer
% AAE 532
%
% load_constants:
%   Load various constants that are often used in orbital analysis.
%   All values rounded to six significant figures.
%%%%%

function [] = load_constants()
    % declare all constants as global variables
    global G R_Sun R_Moon R_Mercury R_Venus R_Earth R_Mars R_Jupiter ...
           R_Saturn R_Uranus R_Neptune R_Pluto mu_Sun mu_Moon ...
           mu_Mercury mu_Venus mu_Earth mu_Mars mu_Jupiter mu_Saturn ...
           mu_Uranus mu_Neptune mu_Pluto r_Moon r_Mercury r_Venus ...
           r_Earth r_Mars r_Jupiter r_Saturn r_Uranus r_Neptune r_Pluto

    G = 6.67408e-20;  % universal gravitational constant, km^3/(kg*s^2)
    
    % mean equatorial radius, km
    R_Sun = 695990;
    R_Moon = 1737.40;
    R_Mercury = 2439.70;
    R_Venus = 6051.80;
    R_Earth = 6378.14;
    R_Mars = 3396.19;
    R_Jupiter = 71492.0;
    R_Saturn = 60268.0;
    R_Uranus = 25559.0;
    R_Neptune = 25559.0;
    R_Pluto = 1195.00;
    
    % gravitational parameter (Gm), km^3/s^2
    mu_Sun = 132712000000;
    mu_Moon = 4902.80;
    mu_Mercury = 22032.1;
    mu_Venus = 324859;
    mu_Earth = 398600;
    mu_Mars = 42828.4;
    mu_Jupiter = 126687000;
    mu_Saturn = 37931300;
    mu_Uranus = 5793970;
    mu_Neptune = 6835110;
    mu_Pluto = 873.767;
    
    % semi-major axis of orbit (km)
    r_Moon = 384400;  % around Earth
    r_Mercury = 57909000;
    r_Venus = 108209000;
    r_Earth = 149650000;
    r_Mars = 227953000;
    r_Jupiter = 779067000;
    r_Saturn = 1426650000;
    r_Uranus = 2873680000;
    r_Neptune = 4492500000;
    r_Pluto = 5868130000;