function TA = mean2trueAnomaly(M,e)
% mean2trueAnomaly computes the true anomaly for a given mean anomaly in  
% the elliptic orbit
% INPUT: n [rad/s] = Mean motion of the orbit
%        e [-]     = eccentricity of the orbit
%        M [rad]   = mean anomaly
% 
% OUTPUT: TA [rad]  = True anomaly of the given time from periapsis

%Author: Zeno Pavanello 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
E  = fzero(@(E) M-E+e*sin(E), 0);         %[rad] eccentric anomaly;
TA = 2*atan2(sqrt(1 + e) * sin(E / 2),sqrt(1 - e) * cos(E / 2));  %[rad] true anomaly;
end

