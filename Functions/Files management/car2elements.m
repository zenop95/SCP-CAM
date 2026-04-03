function x = car2elements(x0, mu)

% -------------------------------------------------------------------------
% INPUTS
% - x0          : initial state in ECI J2000 reference frame [km, km/s]
% - mu
% -------------------------------------------------------------------------
% OUTPUT:
% - x         : vector with orbital elements [km, rad]:
%       - a        : semi-major axis [km]
%       - e        : eccentricity [-]
%       - i        : inclination [rad]
%       - Om       : RAAN (Right Ascension of the Ascending Node) [rad]
%       - om       : periapsis true anomaly [rad]
%       - theta    : true anomaly [rad]
%       - T        : orbital period [s]
%       - M        : mean anomaly [rad]
% -------------------------------------------------------------------------
% Author:   Maria Francesca Palermo, Politecnico di Milano, 27 November 2020
%           e-mail: mariafrancesca.palermo@mail.polimi.it

rr = x0(1:3);
vv = x0(4:6);

r = norm(rr);
v = norm(vv);
k = [0 0 1]';

a = - mu / (v^2 - 2*mu/r);

hh = cross(rr,vv);
h = norm(hh);

ee = cross(vv,hh)/mu - rr/r;
e = norm(ee);

i = acos(hh(3)/h);

N = cross(k,hh) / norm(cross(hh,k));

Om = acos(N(1));
if N(2) < 0
    Om = 2*pi - Om;
end

om = acos(dot(N,ee)/e);
if ee(3) < 0
    om = 2*pi - om;
end

theta = acos(dot(rr,ee)/(r*e));
if dot(rr,vv) < 0
    theta = 2*pi - theta;
end

% period
T = 2*pi*sqrt(a^3/mu);

% eccentric anomaly
E = 2* atan(sqrt( (1-e)/(1+e) )*tan(theta/2));

% mean anomaly
M = E - e*sin(E);

% Output
x(1) = a;
x(2) = e;
x(3) = i;
x(4) = Om;
x(5) = om;
x(6) = theta;
x(7) = T;
x(8) = M;

end

