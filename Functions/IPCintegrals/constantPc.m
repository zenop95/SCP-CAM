function Pc = constantPc(mu,P,R)
% constantIpc finds the Instantaneous Probability of Collision between 
% two space objects considering a constant value of the pdf of The relative
% distance between the two spacecraft, described by a gaussian RV.

% INPUT: mu [m]        = exprected value of the relative position RV in ECI
%        P [m^2]       = covariance matrix of the relative position RV in ECI    
%        R [m]         = combined Hard Body Radius of the two spacecraft
%        R [m]         = combined Hard Body Radius of the two spacecraft
%        e2b           = rotation matrix from EXi to B-space
% 
% OUTPUT: Pc [-]       = Probability of Collision; 
%
% Reference: Alfriend, K., Akella, M., Frisbee, J., Foster, J., Lee, D.-J.,
%            & Wilkins, M. (1999). Probability of Collision Error Analysis. 
%            Space Debris, 1(1), 21–35. 
%            https://doi.org/10.1023/A:1010056509803
%
% Author: Zeno Pavanello, 2022
%--------------------------------------------------------------------------
de        = det(P);
mu        = toColumn(mu);
if abs(de) > 1e-60
    smd = mu'*(P\mu);
    Pc = R^2/(2*sqrt(de))*exp(-smd/2);
    Pc(Pc>1) = 1;
else
    Pc = nan;
    warning(['Impossible to compute the constant pdf PC because the ' ...
        'covariance matrix is not positive definite'])
end

