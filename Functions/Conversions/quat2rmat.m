function R = quat2rmat(q)
% Compute the rotation matrix corresponding to the input quaternion
%
% Usage:
%   R = quat2rmat(q)
%
% Inputs:
%   q: 4xN unit quaternion
% Outputs:
%   R: 3x3xN corresponding rotation matrix

[nrow, ncol] = size(q);
if ~isequal(nrow, 4) && isequal(ncol, 4)
    % warning('Expected input dimensions are 4xN, found %dx%d. Will transpose the input vector.', nrow, ncol);
    q = q';
end
validateattributes(q,{'numeric'},{'nrows',4,'real','finite','nonnan'})
n = size(q, 2);

% Normalize, otherwise the matrix will be wrong
q = quat_normalize(q);

% Reshape the quaternions to be in the 3rd dimension
q = reshape(q,[4 1 size(q,2)]);
% Extract components from input quaternion
w = q(1,1,:); x = q(2,1,:); y = q(3,1,:); z = q(4,1,:);

r = cat(1,...
    1 - 2*(y.^2 + z.^2),    2*(x.*y - w.*z),        2*(x.*z + w.*y),...
    2*(x.*y + w.*z),        1 - 2*(x.^2 + z.^2),    2*(y.*z - w.*x),...
    2*(x.*z - w.*y),        2*(y.*z + w.*x),        1 - 2*(x.^2 + y.^2) );

% Prepare output in the correct format
R = permute(reshape(r, 3, 3, n), [2 1 3]);
end