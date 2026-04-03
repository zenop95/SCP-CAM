function quat_conj = quat_conjugation(q)
% Compute the conjugate of the input quaternion
% quat = [eta; eps] -> quat_conj = [eta; -eps]
%
%  Input:
%   q: 4xN matrix containing M quaternions.
%
%   Output:
%   q*: 4xN matrix containing conjugation of quaternions.

validateattributes(q,{'numeric'},{'nrows',4,'real','finite','nonnan'})

len = size(q,2);
quat_conj   =   reshape(cell2mat(arrayfun(@(i) [q(1,i); -q(2:4,i)],1:len,...
    'UniformOutput',false)'),4,len);
end