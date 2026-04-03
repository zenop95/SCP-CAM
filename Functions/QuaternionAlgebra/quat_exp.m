function [e_q] = quat_exp(q)
%Computes the exponential of a quaternion
% Input:
%   q:       4-by-M matrix containing M quaternions.
%
% Output:
%   e_q:     4-by-M matrix containing exponentials of quaternions.

validateattributes(q,{'numeric'},{'nrows',4,'real','finite','nonnan'})

len = size(q,2);
v_norm = arrayfun(@(i) norm(q(2:4,i)),1:len, 'UniformOutput',true);

% Calculate the exponential
e_q = reshape(cell2mat(arrayfun(@(i) exp(q(1,i))*[cos(v_norm(i));...
    sin(v_norm(i))*q(2:4,i)/v_norm(i)],1:len, 'UniformOutput',false)'),4,len);
% Handle singularities
e_q(isnan(e_q)) = 0;
end

