function qinv = quat_inv(q)
% Compute the inverse of the input quaternion
%
% Input:
%   q:            4-by-M matrix containing M quaternions.
%
% Output:
%   qinv:         4-by-M matrix containing the inverse of quaternions.

validateattributes(q,{'numeric'},{'nrows',4,'real','finite','nonnan'})

len = size(q,2);
qinv  = reshape(cell2mat(arrayfun(@(i)...
    quat_conjugation(q(:,i))./(quat_norm(q(:,i))^2),1:len, 'UniformOutput',false)'),4,len);
end