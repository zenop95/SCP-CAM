function IPC = constantIpc(mu,P,R)
% constantIpc finds the Instantaneous Probability of Collision between 
% two space objects considering a constant value of the pdf of The relative
% distance between the two spacecraft, described by a gaussian RV.

% INPUT: mu [m]        = exprected value of the relative position RV
%        P [m^2]       = covariance matrix of the relative position RV       
%        R [m]         = combined Hard Body Radius of the two spacecraft

% OUTPUT: maxIPC [-] = maximum Instantaneous Probability of Collision; 

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
    IPC = sqrt(2)*R^3/(3*sqrt(pi*de))*exp(-sqrMahalanobis/2);
    IPC(IPC>1) = 1;
else
    IPC = nan;
    warning(['Impossible to compute the constant pdf ICP because the ' ...
        'covariance matrix is not positive definite'])
end

