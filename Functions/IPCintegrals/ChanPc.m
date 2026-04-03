function PoC = ChanPc(r,P,R,varargin)
% computePoC computes the PoC with Chan's method.
%
% INPUT: r     = [-] (2,1)   Relative trajectory 
%        P     = [-] (2,2,1) Positional covariance
%        R     = [-] (1,1)   Hard body radius
%        order = [-] (1,1)   Maximum order of the series
%        
% OUTPUT: PoC = [-] (1,1) Probability of Collision
% 
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
if nargin > 3; order = varargin{1}; else; order = 5; end
sigma_xi   = sqrt(P(1,1));
sigma_zeta = sqrt(P(2,2));
rho        = abs(P(1,2))/(sigma_xi*sigma_zeta);

u          = R^2/(sigma_xi*sigma_zeta*sqrt(1.0-rho^2));
v          = (((r(1)/sigma_xi))^2 + ((r(2)/sigma_zeta))^2 - ...
                2*rho*(r(1)*r(2)/(sigma_xi*sigma_zeta)))/(1-rho^2);
outerSum = 0;

for m = 0:order
    innerSum = 0;
    for k = 0:m
        innerSum = innerSum + u^k/(2^k*factorial(k));
    end
    outerSum = outerSum + v^m/(2^m* ...
                            factorial(m))*(1 - exp(-u/2)*innerSum);
end
PoC = exp(-v/2)*outerSum;

end