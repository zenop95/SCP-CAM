function IPC = cuboidIpc(mu,P,R)
% cuboidIpc computes the Instantaneous Collision Probability (IPC) as the
% integral of the gaussian probability distribution ellipsoid of the 
% relative position of two bodies described by the mean and covariance 
% matrix in input over the Hard Body Sphere of radius R

% INPUT: mu = [m] expected value of the relative position random variable
%        P  = [m^2] Covariance Matrix of the relative position random variable
%        R  = [m] Hard Body Radius (HBR) of the collision

% OUTPUT: Ipc = [-] Computed Instantaneous Collision Probability


% Bibliography: Zhang, S., Fu, T., Chen, D., & Cao, H. (2020). 
%               Satellite instantaneous  collision probability computation 
%               using equivalent volume cuboids. Journal of Guidance, 
%               Control, and Dynamics, 43(9), 1757–1763. 

% Author: Zeno Pavanello, 2022
%--------------------------------------------------------------------------

%Factorize covariance matrix
try L = chol(P)';
    %Scaled gaussian ellipsoid
    muPrime = L\mu;
    A = L'*L;
    %Find orthogonal matrix of eigenvectors of A and its eigenvalues
    [Q,Lambda] = eig(A);
    
    %Rotation Transform and Equivalent cuboid transform (to preserve volume)
    muSecond = Q'*muPrime;
    extremesNew = (pi/6)^(1/3)*R./sqrt(diag(Lambda));
    
    I = nan(3,1);
    %Final integration 
    for i = 1:3
    I(i) = (erf((muSecond(i)+extremesNew(i))/sqrt(2))- ...
            erf((muSecond(i)-extremesNew(i))/sqrt(2)))/2;
    end
    IPC = prod(I);

catch 
    warning('The Covariance Matrix is not positive definite')
    if abs(det(P)) < 1e-14
        warning('The covariance matrix is singular')
    end
    IPC = nan;
end 
end

