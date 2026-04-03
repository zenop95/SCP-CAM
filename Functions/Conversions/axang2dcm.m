function dcm = axang2dcm(axis, ang)
% Compute the direction cosine matrix describing the rotation of the input 
% angle around the input axis according to Rodrigues' rotation formula
%
% Usage:
%   dcm = axang2dcm(axis, angle)
%
% Inputs:
%   axis:  3xN vector defining the rotation axis
%   angle: 1xN rotation angle [rad]
%
% Outputs:
%   dcm:   3x3xN corresponding direction cosine matrix. The provided DCM 
%          is in the pre-multiply format (i.e., v2 = dcm21 * v1).

validateattributes(axis,{'numeric'},{'nrows',3,'real','finite','nonnan'})
validateattributes(ang,{'numeric'},{'nrows',1,'real','finite','nonnan'})

R = axang2rmat(axis, ang);
dcm = R';
end