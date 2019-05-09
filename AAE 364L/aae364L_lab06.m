%%%%%
% Jordan Mayer
% AAE 364L
% Lab 06
%
% Do things
%%%%%

%% Preliminary setup
close all; clear all; format compact; format short; rehash toolbox;
addpath('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 06\Lab 06 Files');
addpath('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 06\Lab 06 Data');
load('lab06');  % load experimental data
svp.time = statevoltage_pitch_benny(:,1);
svp.signals.values = statevoltage_pitch_benny(:,2:3);
statevoltage_pitch_benny = svp;

%% Get performance criteria

[epp_L2e, epp_eMss] = get_errs(errtheta_pitch_benny)
[epy_L2e, epy_eMss] = get_errs(errtheta_yaw_benny)
[spp_L2e, spp_eMss] = get_errs(statepitch_pitch_benny)
[spy_L2e, spy_eMss] = get_errs(statepitch_yaw_benny)
[eyp_L2e, eyp_eMss] = get_errs(errpsi_pitch_benny)
[eyy_L2e, eyy_eMss] = get_errs(errpsi_yaw_benny)
[syp_L2e, syp_eMss] = get_errs(stateyaw_pitch_benny)
[syy_L2e, syy_eMss] = get_errs(stateyaw_yaw_benny)

%% Get plots
close all;

ttl_start = 'Jordan Mayer, AAE 364L, Lab 06, ';
ttl_i = [ttl_start, 'Part (i) Results'];
ttl_ii = [ttl_start, 'Part (ii) Results'];

pitch_lab = 'Pitch Angle, deg';
yaw_lab = 'Yaw Angle, deg';
vp_lab = 'Pitch Voltage, V';
vy_lab = 'Yaw Voltage, V';

create_plot(errtheta_pitch_benny, pitch_lab, ttl_i);
create_plot(errpsi_yaw_benny, yaw_lab, ttl_i);
plot_volt(errvoltage_pitch_benny, [ttl_i, ', Pitch Input']);
plot_volt(errvoltage_yaw_benny, [ttl_i, ', Yaw Input']);

create_plot(statepitch_pitch_benny, pitch_lab, ttl_ii);
create_plot(stateyaw_yaw_benny, yaw_lab, ttl_ii);
plot_volt(statevoltage_pitch_benny, [ttl_ii, ', Pitch Input']);
plot_volt(statevoltage_yaw_benny, [ttl_ii, ', Yaw Input']);

compare_plot(errtheta_pitch_benny, statepitch_pitch_benny, ...
             pitch_lab, [ttl_start, ', Comparison']);
compare_plot(errpsi_yaw_benny, stateyaw_yaw_benny, ...
             yaw_lab, [ttl_start, ', Comparison']);     