function rv = COE2RV(Coe, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	rv = COE2RV(Coe, mu)
% 	
% 	
% 	INPUT:
% 	- coe 		(1x6) orbit element array (L, -, rad x 4)
%     .p       
%     .f       
%     .g      
%     .h    
%     .L    
% 	- mu 		gravitational parameter (L^3/s^2)
% 	
% 	OUTPUT:
% 	- a
% 	- e
% 	- in
% 	- OM
% 	- om
% 	- theta
% 	
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 1
    mu = varargin{1};
else
    mu = 1;
end

elimitpar = 1e-17;

a   = Coe(1);
e   = Coe(2);
i   = Coe(3);
Om  = Coe(4);
om  = Coe(5);
tho = Coe(6);

% Rotation matrix
R = zeros(3,3);
R(1,1) = cos(om)*cos(Om)-sin(om)*cos(i)*sin(Om);
R(2,1) = cos(om)*sin(Om)+sin(om)*cos(i)*cos(Om);
R(3,1) = sin(om)*sin(i);

R(1,2) = -sin(om)*cos(Om)-cos(om)*cos(i)*sin(Om);
R(2,2) = -sin(om)*sin(Om)+cos(om)*cos(i)*cos(Om);
R(3,2) = cos(om)*sin(i);

R(1,3) = sin(i)*sin(Om);
R(2,3) = -sin(i)*cos(Om);
R(3,3) = cos(i);

% In plane Parameters
% ---> CL: 07/10/2010, Matteo Ceriotti: added condition e <= (1+elimitpar)
if e >= (1-elimitpar) && e <= (1+elimitpar) % Parabola
    if nargin < 3
        error('Parabolic case: the semi-latus rectum needs to be provided')
    end
else
    % ---> CL: 07/10/2010, Matteo Ceriotti: Removed if nargin < 3. p shall
    % be computed even if it is given, if it is not a parabola.
    p = a*(1-e^2);     % Value of p in the input is not considered and overwritten with this one
end

r = p/(1+e*cos(tho));
xp = r*cos(tho);
yp = r*sin(tho);
wom_dot = sqrt(mu*p)/r^2;
r_dot   = sqrt(mu/p)*e*sin(tho);
vxp = r_dot*cos(tho)-r*sin(tho)*wom_dot;
vyp = r_dot*sin(tho)+r*cos(tho)*wom_dot;

% 3D cartesian vector
rv = zeros(6,1);
rv(1) = R(1,1)*xp+R(1,2)*yp;
rv(2) = R(2,1)*xp+R(2,2)*yp;
rv(3) = R(3,1)*xp+R(3,2)*yp;

rv(4) = R(1,1)*vxp+R(1,2)*vyp;
rv(5) = R(2,1)*vxp+R(2,2)*vyp;
rv(6) = R(3,1)*vxp+R(3,2)*vyp;

end
