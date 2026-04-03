function dcm = quat2dcm(q)
% Compute the direction cosine matrix corresponding to the input quaternion
%
% Usage:
%   dcm = quat2rmat(q)
%
% Inputs:
%   q: 4xN unit quaternions
% Outputs:
%   dcm: 3x3xN corresponding direction cosine matrix

% Rotating input to the quat2rmat function
[nrow, ncol] = size(q);
if ~isequal(nrow, 4) && isequal(ncol, 4)
    % warning('Expected input dimensions are 4xN, found %dx%d. Will transpose the input vector.', nrow, ncol);
    q = q';
end
validateattributes(q,{'numeric'},{'nrows',4,'real','finite','nonnan'})

% Output a transposed rotation matrix/matrices array
dcm = permute(quat2rmat(q), [2 1 3]);
end