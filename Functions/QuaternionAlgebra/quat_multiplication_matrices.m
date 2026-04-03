function [matM,matN] = quat_multiplication_matrices(quat_1,quat_2)
% Compute the matrices involved in the multiplication (composition) of the two input quaternions
% quat_1 = [eta_1; eps_1], quat_2 = [eta_2; eps_2] -> M(quat_1) and N(quat_2) sunch that
% quat_res = quat_multiplication(quat_1,quat_2) = M(quat_1)*quat_2 = N(quat_2)*quat_1
%
% INPUT
%   quat_1       [4x1]   quaternion (possibly as column vector)
%   quat_2       [4x1]   quaternion (possibly as column vector)
% OUTPUT
%   matM         [4x3]   matrix M involved in quaternions multiplication
%   matN         [4x3]   matrix N involved in quaternions multiplication

quat_1      =   quat_reshape(quat_1); 
quat_2      =   quat_reshape(quat_2); 

matM        =   zeros(4,4);
matN        =   zeros(4,4);
skew_1      =   [0 -quat_1(4) quat_1(3); quat_1(4) 0 -quat_1(2); -quat_1(3) quat_1(2) 0];
skew_2      =   [0 -quat_2(4) quat_2(3); quat_2(4) 0 -quat_2(2); -quat_2(3) quat_2(2) 0];
matM        =   [quat_1(1) -quat_1(2:4)'; quat_1(2:4) quat_1(1)*eye(3)+skew_1];
matN        =   [quat_2(1) -quat_2(2:4)'; quat_2(2:4) quat_2(1)*eye(3)-skew_2];
end