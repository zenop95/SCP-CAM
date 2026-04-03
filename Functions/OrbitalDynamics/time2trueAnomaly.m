function TA = time2trueAnomaly(n,e,t)
% time2trueAnomaly computes the true anomaly for a given time in the 
% elliptic orbit
% INPUT: n [rad/s] = Mean motion of the orbit
%        e [-]     = eccentricity of the orbit
%        t [s]     = time from perigee
% 
% OUTPUT: TA [rad]  = True anomaly of the given time from periapsis

%Author: Zeno Pavanello 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
M  = n*t;                      %[rad] mean anomaly;
TA = mean2trueAnomaly(n,e,M,t);  %[rad] true anomaly;
TA = mod(TA,2*pi);
end

