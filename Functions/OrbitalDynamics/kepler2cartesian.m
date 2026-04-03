function x = kepler2cartesian(a,e,RA,inc,w,TA,varargin)
% This function computes the state vector (r,v) expressed in ECI
% from the classical orbital elements (coe).
%
% Usage: [r, v] = kepler2cartesian(a,e,RA,inc,w,TA, mu)
%
% Inputs
% coe - structure containing the orbital elements [a e RA incl w TA]
%    where
%           a    = semi major axis (km\m)
%           e    = eccentricity
%           RA   = right ascension of the ascending node (rad)
%           incl = inclination of the orbit (rad)
%           w    = argument of perigee (rad)
%           TA   = true anomaly (rad)
% Optional inputs:
%   mu 	- [kmˆ3; sˆ2] gravitational parameter. Defaults to 398600.4418.
%   
% Outputs
%   r	- [km] position vector in the geocentric equatorial frame
%   v 	- [km/s] velocity vector in the geocentric equatorial frame

% Used variables:
% R3_w - Rotation matrix about the z-axis through the angle w
% R1_i - Rotation matrix about the x-axis through the angle i
% R3_W - Rotation matrix about the z-axis through the angle RA
% Q_pX - Matrix of the transformation from perifocal to
%        geocentric equatorial frame
% rp   - position vector in the perifocal frame (km)
% vp   - velocity vector in the perifocal frame (km/s)
%
% Equation numbers refer to Curtis, ed. 3.

% Authors: Davide Vertuani, Zeno Pavanello 2021
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

if nargin == 7 
    mu = varargin{1}; 
else 
    mu = 3.986004418e14; 
end

h = sqrt(a*mu*(1-e^2));

% Equations 4.37 and 4.38 (rp and vp are column vectors):
rp = (h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] ...
    + sin(TA)*[0;1;0]);
vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);

% Equation 4.39:
R3_W = [ cos(RA)  sin(RA)  0
    -sin(RA)  cos(RA)  0
    0        0     1];

% Equation 4.40:
R1_i = [1       0      0
        0   cos(inc)  sin(inc)
        0  -sin(inc)  cos(inc)];

% Equation 4.41:
R3_w = [ cos(w)  sin(w)  0
        -sin(w)  cos(w)  0
         0       0     1];

% Equation 4.44:
Q_pX = R3_W'*R1_i'*R3_w';

% Equations 4.46 (r and v are column):
r = Q_pX*rp;
v = Q_pX*vp;
x = [r; v];
end