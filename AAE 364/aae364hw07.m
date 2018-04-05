%%%%%
% AAE 364
% HW 07
% Jordan Mayer
%
% For various control systems, verify the root loci plotted by hand in part
% a). That is, for each problem:
%   1. Generate L(s) using MATLAB's tf([num], [den]) function.
%   2. Create the root locus plot in a new figure with figure(num) and
%   rlocus(L).
%   3. Give the plot a distinctive title and scale axes for comparison.
%%%%%

close all;  % close all open figures

%%% Problem 1: B-6-1 %%%
L = tf([0 1 1], [1 0 0]);
figure(1); rlocus(L);
title('Root Locus for Problem 1: B-6-1');
axis([-3 3 -3 3]);pbaspect([1 1 1]);

%%% Problem 1: B-6-2 %%%
L = tf([0 0 0 0 1], [1 5 9 5 0]);
figure(2); rlocus(L);
title('Root Locus for Problem 1: B-6-2');
axis([-3 3 -3 3]);pbaspect([1 1 1]);

%%% Problem 1: B-6-6 %%%
L = tf([0 0 1 9], [1 4 11 0]);
figure(3); rlocus(L);
title('Root Locus for Problem 1: B-6-6');
axis([-9 9 -9 9]);pbaspect([1 1 1]);

%%% Problem 1: B-6-8 %%%
L = tf([0 0 0 1], [1 4 8 0]);
figure(4); rlocus(L);
title('Root Locus for Problem 1: B-6-8');
axis([-5 5 -5 5]);pbaspect([1 1 1]);

%%% Problem 2 %%%
L = tf([0 0 0 0 1], [1 4 11 14 10]);
figure(5); rlocus(L);
title('Root Locus for Problem 2');
axis([-3 3 -3 3]);pbaspect([1 1 1]);

%%% Problem 3 %%%
L = tf([0 0 0 2 2], [1 7 10 0 0]);
figure(6); rlocus(L);
title('Root Locus for Problem 3');
axis([-7 7 -7 7]);pbaspect([1 1 1]);