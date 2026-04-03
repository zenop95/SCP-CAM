function [x] = propKep(x0,t,varargin)
% PROPKEP Uses the true anomaly formulation to propagate an unperturbed
% keplerian orbit from state x0 to the final state x(t(end)).
%
% INPUT: x0   = [m][m/s] 6x1 initial state of the satellite
%        t    = [s] Nx1 vector of the time nodes
% 
% OUTPUT: x   = [m][m/s] 6xN matrix of the state of the orbit at each node
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
if nargin == 3 
    mu = varargin{1}; 
else 
    mu = 3.986004418e14; 
end
N      = length(t);
coe    = cartesian2kepler(x0);
theta  = nan(N,1);
tt0    = trueAnomaly2time(coe.n,coe.ecc,coe.theta);
x      = nan(6,N);
x(:,1) = x0;
% Propagate using true anomaly
for i = 2:N
    theta(i) = time2trueAnomaly(coe.n,coe.ecc,tt0+t(i));
    x(:,i)   = kepler2cartesian(coe.a,coe.ecc,coe.RAAN, ...
                                coe.inc,coe.w,theta(i),mu);
end
end

