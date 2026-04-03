function vout = quatrot(v, q)
%  Rotate a vector of the rotation described by a quaternion.
%
% Usage:
%   vout = quatrot(v, q)
%
% Input:
%   v: 3xM vector to be rotated.
%   q: 4xM matrix containing M quaternions.
%
% Output:
%   vout: 3xM matrix containing M rotated vectors.

validateattributes(v,{'numeric'},{'nrows',3,'real','finite','nonnan'})
validateattributes(q,{'numeric'},{'nrows',4,'real','finite','nonnan'})

% Check number of elements
n = size(q,2);
assert(isequal(size(v,2), n), ...
    'Input matrices v and q must have the same number of elements.');

q = quat_normalize(q);
vq = [zeros(1,n); v];

% Controintuitive, as the order is shifted. Must be checked.
vr = cell2mat(arrayfun(@(k) ...
    quat_multiplication(quat_conjugation(q(:,k)), quat_multiplication(vq(:,k), q(:,k))),...
    1:n, 'UniformOutput', false));

vout = vr(2:4,:);
end