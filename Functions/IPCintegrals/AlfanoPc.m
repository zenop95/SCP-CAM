function PoC = AlfanoPc(r,P,HBR)
% AlfanoPc computes the PoC with Alfano's method.
%
% INPUT: r = [-] (3,1)   Relative trajectory 
%        P = [-] (3,3,1) Positional covariance
%        R = [-] (1,1)   Hard body radius
%        
% OUTPUT: PoC = [-] (1,1) Probability of Collision
% 
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
sigma_x = sqrt(P(1,1));
sigma_z = sqrt(P(2,2));
rho = P(1,2)/(sigma_x*sigma_z);
theta = 1/2*atan(2*rho*sigma_x*sigma_z/(sigma_x*sigma_x-sigma_z*sigma_z));

if (sigma_z > sigma_x)
    theta = theta + pi/2;
end

R = [cos(theta) sin(theta); -sin(theta) cos(theta)];

C  = R*P*R';
xx = R*r;
xm = xx(1);
zm = xx(2);

sigmax = sqrt(C(1,1));
sigmaz = sqrt(C(2,2));

% n = floor(5*HBR/min(sqrt(sigmaz),vnorm(r)));
n   = 30;
PoC = 0;

for i = 1:n
    aux1 = ( zm + 2.0*HBR/n*sqrt((n-i)*i) )/(sigmaz*sqrt(2.0));
    aux2 = (-zm + 2.0*HBR/n*sqrt((n-i)*i) )/(sigmaz*sqrt(2.0));
    aux3 = -(HBR*(2.0*i-n)/n + xm)^2 /(2*sigmax*sigmax);
    PoC  = PoC + (erf(aux1)+erf(aux2))*exp(aux3);
end

PoC = HBR*2/(sqrt(8.0*atan(1.0)*4)*sigmax*n)*PoC;

end