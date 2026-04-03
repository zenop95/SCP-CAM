function R = axang2rmat(axis, ang)
% Compute the rotation matrix describing the rotation of the input angle
% around the input axis according to Rodrigues' rotation formula
%
% Usage:
%   R = axang2rmat(axis, angle)
%
% Inputs:
%   axis:  3x1 vector defining the rotation axis
%   angle: 1x1 rotation angle [rad]
% Outputs:
%   R:     3x3 corresponding rotation matrix

quat = axang2quat(axis, ang);
R    = quat2rmat(quat);
end