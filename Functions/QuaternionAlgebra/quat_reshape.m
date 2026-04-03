function q = quat_reshape(q)
% Reshape the input quaternion as a column vector
% 
% INPUT
%   quat: 4xN or Nx4 quaternions array
% OUTPUT
%   quat_res: 4xN quaternion shaped as column vector

[nrow, ncol] = size(q);
assert(nrow~=ncol, 'Ambiguous case, input matrix is 4x4, cannot reshape.')

if size(q, 1) == 4; return; end
q = q';
end