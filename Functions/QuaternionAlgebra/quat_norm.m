function norm = quat_norm(q)
% Compute the norm of the input quaternion
%
% INPUT
%   q:            4-by-M matrix containing M quaternions.
% OUTPUT
%   norm:         1-by-M array containing M quaternion norms.

validateattributes(q,{'numeric'},{'nrows',4,'real','finite','nonnan'})

len = size(q,2); 
norm = reshape(cell2mat(arrayfun(@(i) sqrt(dot(q(:,i),q(:,i))), ...
    1:len, 'UniformOutput',false)'),1,len);
end