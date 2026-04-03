function nsCoe = nonSingularCoe(coe, varargin)
% coeFromCartesian converts the ECI position and velocity vector of the spacecraft
% into the definition of the orbit in terms of six Keplerian parameters
%
% Inputs:
%   coe: structure containing the orbital elements [a e RA inc w TA]
%    where
%           a [m]     = semi major axis 
%           e [-]      = eccentricity
%           RA [rad]   = right ascension of the ascending node 
%           incl [rad] = inclination of the orbit 
%           w [rad]    = argument of perigee
%           TA [rad]   = true anomaly
%           n [rad/s]  = mean motion
% Outputs:
% nsCoe: structure containing the non singular orbital elements 
%    where
%           eVec  [] = eccentricity vector
%           iVec  [] = inclination vector
%           drift [m; []] = drift vector 


% Author: Zeno Pavanello, 2021
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

a     = coe.a;
ecc   = coe.ecc;
inc   = coe.inc;
RAAN  = coe.RAAN;
w     = coe.w;
theta = coe.theta;
E = 2*atan(sqrt((1-ecc)/(1+ecc))*tan(theta/2));  %[rad] eccentric anomaly;
M = E-ecc*sin(E);                           %[rad] mean anomaly;

eVec  = [ecc*cos(RAAN+w); ecc*sin(RAAN+w)];
iVec  = [inc*cos(RAAN); inc*sin(RAAN)]; 
drift = [a; RAAN + w + M - theta];

nsCoe = struct(...
           'eVec',     eVec, ...
           'iVec',   iVec, ...
           'drift',  drift ...
            );
end