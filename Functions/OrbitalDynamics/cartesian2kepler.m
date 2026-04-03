function coe = cartesian2kepler(x, varargin)
% coeFromCartesian converts the ECI position and velocity vector of the spacecraft
% into the definition of the orbit in terms of six Keplerian parameters
%
% Inputs:
%   r [m] : orbital position in ECI reference frame
%   v [m/s] : orbital velocity in ECI reference frame
% Optional inputs:
%   mu 	- [m^3 sˆ2] gravitational parameter. Defaults to 3.986004418e14.
%
% Outputs:
% coe - structure containing the orbital elements [a e RA incl w TA]
%    where
%           a [m]     = semi major axis 
%           e [-]      = eccentricity
%           RA [rad]   = right ascension of the ascending node 
%           incl [rad] = inclination of the orbit 
%           w [rad]    = argument of perigee
%           TA [rad]   = true anomaly
%           n [rad/s]  = mean motion

% Authors: Davide Vertuani, Zeno Pavanello, 2021
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
r = x(1:3);
v = x(4:6);
r = reshape(r, 3, []);
v = reshape(v, 3, []);
validateattributes(r, {'numeric'}, {'numel', 3});
validateattributes(v, {'numeric'}, {'numel', 3});

mu = 398600.4418;
if nargin > 1; mu = varargin{1}; end

r_n     = norm(r);              % norm of Radial distance
v_n     = norm(v);              % norm of Speed
vr      = dot(r,v)/r_n;         % Radial velocity
H       = cross(r,v);           % Specific angular momentum
h       = norm(H);              % Magnitude of the specific angular momentum
inc     = acos(H(3)/h);         % Inclination

if isnan(inc); inc = 0; end

N   = cross([0;0;1],H);      	% Node line vector
Nn   = norm(N);                  % Magnitude of N

 % Right ascension of the ascending node
if(N(2) >= 0); RA = acos(N(1)/Nn);
else; RA = 2*pi - acos(N(1)/Nn);
end

if isnan(RA); RA = 0; end

ev = 1/mu*((v_n^2-mu/r_n)*r-r_n*vr*v);      % Eccentricity vector
e  = norm(ev);                              % Eccentricity
a = h^2/(mu*(1-e^2));                       % Semi-major axis

% Argument of perigee,
if(ev(3) >= 0); w = acos(dot(N,ev)/(Nn*e));
else; w = 2*pi - acos(dot(N,ev)/(Nn*e));
end

if isnan(w); w = 0; end

% True anomaly
if(vr >= 0); TA = acos(dot(ev,r)/(r_n*e));
else; TA = 2*pi - acos(dot(ev,r)/(r_n*e));
end

if isnan(TA); TA = 0; end

% Mean motion
n = sqrt(mu/a^3);

coe = struct(...
           'a',     a, ...
           'ecc',   e, ...
           'RAAN',  RA, ...
           'inc',   inc, ...
           'w',     w, ...
           'theta', TA, ...
           'n',     n ...
            );
end