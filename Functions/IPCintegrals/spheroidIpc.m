function IPC = spheroidIpc(mu,P,R)
% spheroidIpc computes the Instantaneous Collision Probability (IPC) as the
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
    
    %Rotation Transform and Equivalent spheroid transform (to preserve volume)
    muSecond = Q'*muPrime;
    R_new    = (prod(R./sqrt(diag(Lambda))))^(1/3);
    
    %Final integration 
    funHndl = @(x,y,z) 1/(2*pi)^(3/2)*exp(-0.5* ...
                ((x-muSecond(1)).^2+(y-muSecond(2)).^2+(z-muSecond(3)).^2));

    xmin    =         -R_new;
    xmax    =          R_new;
    ymin    = @(x)    -sqrt(R_new^2 - x.^2);
    ymax    = @(x)     sqrt(R_new^2 - x.^2);
    zmin    = @(x,y)  -sqrt(R_new^2 - x.^2 - y.^2);
    zmax    = @(x,y)   sqrt(R_new^2 - x.^2 - y.^2);

    IPC  = integral3(funHndl,xmin,xmax,ymin,ymax,zmin,zmax, ...
                    "AbsTol",1e-10,"RelTol",1e-6);
catch 
    warning('The Covariance Matrix is not positive definite')
    if abs(det(P)) < 1e-14
        warning('The covariance matrix is singular')
    end
    IPC = nan;
end 
end

