function t = trueAnomaly2time(n,e,theta,varargin)
% trueAnomaly2time computes the time from the perigee for a given true 
% anomaly in the elliptic orbit.
% Computes only positive times included between 0 and T where T is the
% orbital period
% INPUT: n [rad/s] = Mean motion of the orbit
%        e [-]     = eccentricity of the orbit
%        theta [rad]  = True anomaly
% 
% OUTPUT: t [s] = time from periapsis of the given true anomaly

% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
t0 = 0;                                   %[s] offset time
T = 2*pi/n;                               %[rad/s] orbital period
if nargin > 3 
    theta0 = varargin{1};
    E0 = 2*atan(sqrt((1-e)/(1+e))*tan(theta0/2));  %[rad] eccentric anomaly;
    M0 = E0-e*sin(E0);                             %[rad] mean anomaly;
    t0 = M0/n;                                     %[s] time
end

E = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2));  %[rad] eccentric anomaly;
M = E-e*sin(E);                              %[rad] mean anomaly;
t = M/n - t0;                                %[s] time
% Eliminate negative times and times over one period
t = mod(t,T);
end

