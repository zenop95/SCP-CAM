function x = propKepOde(x0,u,t,varargin)
% PROPKEPODE Uses solves the ode to propagate the
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
if nargin == 4 
    mu = varargin{1}; 
else 
    mu = 3.986004418e14; 
end
options = odeset('AbsTol',1e-11,'RelTol',1e-12);
x0      = [x0; u];
str     = ode78(@(t,y) grav(y,mu),[0,t],x0,options);
x       = str.y(1:6,end);

end

function dy = grav(y,muS)
    R = y(1:3);
    V = y(4:6);
    uRtn = y(7:9);
    dcm = rtn2eci(R,V);
    u = dcm*uRtn;
    Rmag = norm(R);
    rDot = V;
    vDot = - muS * R / Rmag^3;
    dy = [rDot; vDot + u; zeros(3,1)];
end