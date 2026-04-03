function q = dcm2quat(dcm)
% Compute the quaternion corresponding to the input direction cosine matrix
%
% Usage:
%   q = dcm2quat(dcm)
%
% Inputs:
%   R: 3x3xN direction cosine matrix, in the form to be premultiplied.
%
% Outputs:
%   q: 4xN corresponding  unit quaternion

validateattributes(dcm,{'numeric'},{'size',[3,3,NaN], 'real','finite','nonnan'})

% Transpose each dcm to obtain rotation matrices
R = permute(dcm, [2 1 3]);

% Compute the quaternions
q = rmat2quat(R);
end