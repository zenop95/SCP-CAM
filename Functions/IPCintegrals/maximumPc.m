function Pc = maximumPc(mu,P,R)
% maximumPc finds the maximum Probability of Collision between 
% two space objects. The relative distance between the two spacecraft is
% described by a gaussian RV.

% INPUT: mu [m]        = exprected value of the relative position RV
%        P [m^2]       = covariance matrix of the relative position RV       
%        R [m]         = combined Hard Body Radius of the two spacecraft

% OUTPUT: IPC [-] = maximum Instantaneous Probability of Collision; 

% Reference: Alfriend, K., Akella, M., Frisbee, J., Foster, J., Lee, D.-J.,
%            & Wilkins, M. (1999). Probability of Collision Error Analysis. 
%            Space Debris, 1(1), 21–35. 
%            https://doi.org/10.1023/A:1010056509803

% Author: Zeno Pavanello, 2022
%--------------------------------------------------------------------------
de = det(P);
if abs(de) > 1e-70
    mu = toColumn(mu);
    sqrMahalanobis = mu'*(P\mu);
    Pc = R^2/(exp(1)*sqrt(de)*sqrMahalanobis);
    Pc(Pc>1) = 1;
else
    Pc = nan;
    warning(['Impossible to compute the maximum Pc because the ' ...
        'covariance matrix is not positive definite'])
end

