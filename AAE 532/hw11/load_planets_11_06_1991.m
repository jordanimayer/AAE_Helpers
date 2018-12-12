%%%%%
% Jordan Mayer
% AAE 532
%
% load_planets_11_06_991:
%   Load planetary characteristics for Earth, Jupiter, Saturn, and
%   Neptune for 11/06/1991.
%%%%%

function [] = load_planets_11_06_1991()
    
    % declare as global variables
    global a_Earth a_Jupiter a_Saturn a_Neptune ...
           e_Earth e_Jupiter e_Saturn e_Neptune ...
           i_Earth i_Jupiter i_Saturn i_Neptune ...
           raan_Earth raan_Jupiter raan_Saturn raan_Neptune ...
           aop_Earth aop_Jupiter aop_Saturn aop_Neptune ...
           M0_Earth M0_Jupiter M0_Saturn M0_Neptune;

    % semi-major axis, km
    a_Earth = 149597927.0;
    a_Jupiter = 778328370.0;
    a_Saturn = 1426990810.0;
    a_Neptune = 4496638040.0;
    
    % eccentricity
    e_Earth = .016712542;
    e_Jupiter = .048487326;
    e_Saturn = .055571251;
    e_Neptune = .008599220;
    
    % inclination, rad
    i_Earth = deg2rad(.001051926);
    i_Jupiter = deg2rad(1.30337171);
    i_Saturn = deg2rad(2.48804735);
    i_Neptune = deg2rad(1.77010817);
    
    % right ascension from ascending node, rad
    raan_Earth = deg2rad(0.0);
    raan_Jupiter = deg2rad(100.435428);
    raan_Saturn = deg2rad(113.677893);
    raan_Neptune = deg2rad(131.780238);
    
    % argument of periapsis, rad
    aop_Earth = deg2rad(102.914377);
    aop_Jupiter = deg2rad(-86.2698735);
    aop_Saturn = deg2rad(-20.6752351);
    aop_Neptune = deg2rad(-87.1124809);
    
    % mean anomaly, rad
    M0_Earth = deg2rad(-58.0848626);
    M0_Jupiter = deg2rad(132.596142);
    M0_Saturn = deg2rad(-142.586443);
    M0_Neptune = deg2rad(-117.599104);
end
