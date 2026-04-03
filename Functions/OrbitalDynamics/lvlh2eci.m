function [DCM_IO, w_IO] = lvlh2eci(r, v)
% Calculates the cosine direction matrix to rotate the axis from the
% initial orbital reference frame LVLH to ECI.
% Usage:
%   dcm = eci2lvlh(r, v)
% Inputs:
%   - r: the orbital radius in ECI coordinates
%   - v: the spacecraft velocity in ECI coordinates
%   Note that r and v may be expressed in any measurement unit, as they
%   will be converted to versors; only their orientation is relevant.
% Outputs:
%   - dcm: a 3x3 matrix rotating LVLH into ECI; it might be employed as
%   follows: 
%           v_i = dcm_io * v_o
%   with v_i being a vector in the ECI frame, v_o in the orbital LVLH one.

%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------


% Ensure proper orientation; the DCM might be transposed otherwise
r = reshape(r, [1, 3]);
v = reshape(v, [1, 3]);

% Z defined along the orbital radius with opposite direction
z1 = -r; 
z  = normalize(z1, 'norm');

% Y defined perpendicular to orbital speed and position vector
y1 = cross(v,r);
y  = normalize(y1,'norm');

% X defined perpendicular to Y and Z (i.e., in the direction of speed).
x1 = cross(y, z);         
x  = normalize(x1, 'norm');

%DCM rotating Orbital frame into ECI frame
DCM_IO = [x; y; z]';

if abs(det(DCM_IO)-1) > 1e-8
    warning('DCM is bad conditioned')
end

r = reshape(r, [3,1]);
v = reshape(v, [3,1]);

%angular rate of Orbital wrt to ECI expressed in ECI
w_OI = cross(r,v)/norm(r)^2;

%angular rate of ECI wrt to Orbital expressed in Orbital frame
w_IO = -DCM_IO'*w_OI;
end