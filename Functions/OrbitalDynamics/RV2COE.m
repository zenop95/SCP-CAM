function Coe = RV2COE(rv, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Coe = RV2COE(rv, mu)
% 	
% 	
% 	INPUT:
% 	- a
% 	- e
% 	- in
% 	- OM
% 	- om
% 	- theta
% 	- mu 		gravitational parameter (L^3/s^2)
% 	
% 	OUTPUT:

% 	- coe 		(1x6) orbit element array (L, -, rad x 4)
%     .p       
%     .f       
%     .g      
%     .h    
%     .L    
% 	
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    mu = varargin{1};
else
    mu = 1;
end

elimit = 0.00000001;

r = rv(1:3);
r = r';
v = rv(4:6);
v = v';

nr = sqrt(r(1)^2+r(2)^2+r(3)^2);    % Norm of r

% Angular momentum vector: h = cross(r,v) 
h  = [r(2)*v(3)-r(3)*v(2),r(3)*v(1)-r(1)*v(3),r(1)*v(2)-r(2)*v(1)]; 
nh = sqrt(h(1)^2+h(2)^2+h(3)^2);    % Norm of h

% Inclination
i = acos(h(3)/nh);

% Line of nodes vector
if i ~= 0 && i ~= pi    % Orbit is out of xy plane
    
    % n = cross([0 0 1],h);
    % Normalisation of nv to 1: n = n/norm(n)
    %                           n = n/sqrt(n(1)^2+n(2)^2+n(3)^2);
    n = [-h(2),h(1),0]/sqrt(h(1)^2+h(2)^2);

else                    % Orbit is in xy plane
    
    % Arbitrary choice of n
    n = [1,0,0]; 
%     warning('spaceToolbox:car2kep:planarOrbit','Planar orbit. Arbitrary choice of Omega = 0.');
end

% Argument of the ascending node
Om = acos(n(1));
if n(2) < 0
    Om = mod(2*pi-Om,2*pi);
end

% Parameter
p = nh^2/mu;

% Eccentricity vector: ev = 1/mu*cross(v,h) - r/nr
% ev  = 1/mu*[v(2)*h(3)-v(3)*h(2),v(3)*h(1)-v(1)*h(3),v(1)*h(2)-v(2)*h(1)] - r/nr; 
ev = 1/mu*cross(v,h) - r/nr;
e = sqrt(ev(1)^2+ev(2)^2+ev(3)^2);    % Eccentricity (norm of eccentricity vector)

% Argument of the pericentre
if e<elimit     % Circular orbit
    
    % Arbitrary eccentricity vector
    ev = n;
    ne = 1; % ne = norm(ev)

else                    % Non circular orbit
    ne = e;
end         

om = acos(min(max((n(1)*ev(1) + n(2)*ev(2) + n(3)*ev(3)) / ne,-1),1)); % acos(dot(n,ev)/ne)

if dot(h,cross(n,ev)) < 0 
    om = mod(2*pi-om,2*pi);
end

% Semi-major axis
a = p/(1-e^2);

% True anomaly: acos(dot(ev,r)/ne/nr);
th = acos(min(max((ev(1)*r(1)+ev(2)*r(2)+ev(3)*r(3))/ne/nr,-1),1));

if dot(h,cross(ev,r)) < 0
    % the condition dot(r,cross(h,ev)) < 0 works in the same way
    th = mod(2*pi-th,2*pi);
end

Coe = [a,e,i,Om,om,th]';