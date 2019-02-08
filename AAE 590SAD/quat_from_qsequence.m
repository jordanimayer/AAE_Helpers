%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% quat_from_qsequence:
%   Calculate quaternion from frame A to frame B from quaternions from
%   A to intermediate frame B' and B' to B. Note that vector components of
%   all quaternions will be in B' frame (including final quaternion, which
%   must then be transformed back to A-frame/B-frame).
%
% Inputs:
%   A_q_Bprime: quaternion from A to B', vector in B' frame
%   Bprime_q_B: quaternion from B' to B, vector in B' frame
%
% Outputs:
%   A_q_B: quaternion from A to B, vector in B' frame
%%%%%

function [A_q_B] = quat_from_qsequence(A_q_Bprime, Bprime_q_B)
    A_q_Bprime_vec = A_q_Bprime(1:3);
    A_q_Bprime_4 = A_q_Bprime(4);
    Bprime_q_B_vec = Bprime_q_B(1:3);
    Bprime_q_B_4 = Bprime_q_B(4);
    
    A_q_B_vec = A_q_Bprime_vec * Bprime_q_B_4 + ...
                Bprime_q_B_vec * A_q_Bprime_4 + ...
                cross(Bprime_q_B_vec, A_q_Bprime_vec);
    A_q_B_4 = A_q_Bprime_4 * Bprime_q_B_4 - ...
              dot(A_q_Bprime_vec, Bprime_q_B_vec);
          
    A_q_B = [A_q_B_vec, A_q_B_4];
end