function quat = axang2quat(axis, ang)
% Compute the unit quaternion describing the rotation
% of the input angle around the input axis
%
% Usage:
%   quat = axang2quat(axis, angle)
%
% Inputs:
%   axis:  3xN vector defining the rotation axis
%   angle: 1xN rotation angles [rad]
% OUTPUT
%   quat:  4xN corresponding unit quaternions

validateattributes(axis,{'numeric'},{'nrows',3,'real','finite','nonnan'})
validateattributes(ang,{'numeric'},{'nrows',1,'real','finite','nonnan'})

if ~isequal(size(axis, 2), size(ang, 2))
   error('''axis'' and ''angle'' sizes do not match.'); 
end

quat = [cos(ang/2); axis.*sin(ang/2)];
end