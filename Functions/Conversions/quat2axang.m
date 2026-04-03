function [axis, angle] = quat2axang(q)
% Compute the rotation axis and angle corresponding to the input quaternion
%
% Usage:
%   [axis, angle] = quat2axang(q)
%
% Inputs:
%   q: 4x1 unit quaternion
%
% Outputs:
%   axis:  3x1 vector defining the rotation axis
%   angle: 1x1 rotation angle [rad]

q = q/norm(q);
angle = 2*acos(q(1));

if angle == 0
    axis = [1;0;0];
else
    axis = [q(2,:); q(3,:); q(4,:)] ./ sin(angle/2);
end
end