%%%%%
% Jordan Mayer
% AAE 564
% HW 04
%
% Obtain RREF of matrix to solve system of linear equations
%%%%%

Ab = [1, -1, 2, 1, 5;
      1, -1, 1, 0, 3;
      -2, 2, 0, 2, -2;
      2, -2, -1, -3, 0];
Ab_hat = rref(Ab);

fprintf('rref([A b]) = \n');
disp(Ab_hat);