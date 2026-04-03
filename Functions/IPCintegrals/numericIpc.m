function IPC = numericIpc(mu,P,R)
% numericIcp computes the Instantaneous Collision Probability (IPC) as the
% integral of the gaussian probability distribution ellipsoid of the 
% relative position of two bodies described by the mean and covariance 
% matrix in input over the Hard Body Sphere of radius R

% INPUT: mu = [m] expected value of the relative position random variable
%        P  = [m.^2] Covariance Matrix of the relative position random variable
%        R  = [m] Hard Body Radius (HBR) of the collision

% OUTPUT: 
%        Ipc = [-] Computed Instantaneous Collision Probability


% Bibliography: Núñez Garzón, U. E., & Lightsey, E. G. (2022). Relating 
%               Collision Probability and Separation Indicators in 
%               Spacecraft Formation Collision Risk Analysis. Journal of 
%               Guidance, Control, and Dynamics, 1–16. 
%               https:././doi.org./10.2514./1.g005744

% Author: Zeno Pavanello, 2022
%--------------------------------------------------------------------------
validateattributes(mu,{'double'},{'vector','numel',3})
validateattributes(P,{'double'},{'2d','nrows',3,'ncols',3})
validateattributes(R,{'double'},{'scalar','positive'})
mu      = toColumn(mu);

C       = 1/sqrt((2*pi).^3*det(P));
mu1     = mu(1);
mu2     = mu(2);
mu3     = mu(3);
P1      = P(1,1);
P2      = P(2,2);
P3      = P(3,3);
P12     = P(1,2);
P13     = P(1,3);
P23     = P(2,3);

funHndl = @(x,y,z)  C.*exp(-0.5*((mu3 - z).*(((P12.^2 - P1.*P2).* ...
    (mu3 - z))./(P3.*P12.^2 - 2.*P12.*P13.*P23 + P2.*P13.^2 + P1.*P23.^2 - ...
    P1.*P2.*P3) + ((P2.*P13 - P12.*P23).*(mu1 - x))./(P3.*P12.^2 - 2.* ...
    P12.*P13.*P23 + P2.*P13.^2 + P1.*P23.^2 - P1.*P2.*P3) + ((P1.*P23 - ...
    P12.*P13).*(mu2 - y))./(P3.*P12.^2 - 2.*P12.*P13.*P23 + P2.*P13.^2 + ...
    P1.*P23.^2 - P1.*P2.*P3)) + (mu2 - y).*(((P13.^2 - P1.*P3).*(mu2 - ...
    y))./(P3.*P12.^2 - 2.*P12.*P13.*P23 + P2.*P13.^2 + P1.*P23.^2 - P1.* ...
    P2.*P3) + ((P3.*P12 - P13.*P23).*(mu1 - x))./(P3.*P12.^2 - 2.*P12.* ...
    P13.*P23 + P2.*P13.^2 + P1.*P23.^2 - P1.*P2.*P3) + ((P1.*P23 - P12.* ...
    P13).*(mu3 - z))./(P3.*P12.^2 - 2.*P12.*P13.*P23 + P2.*P13.^2 + P1.* ...
    P23.^2 - P1.*P2.*P3)) + (mu1 - x).*(((P23.^2 - P2.*P3).*(mu1 - x))./ ...
    (P3.*P12.^2 - 2.*P12.*P13.*P23 + P2.*P13.^2 + P1.*P23.^2 - P1.*P2.*P3) + ...
    ((P3.*P12 - P13.*P23).*(mu2 - y))./(P3.*P12.^2 - 2.*P12.*P13.*P23 + ...
    P2.*P13.^2 + P1.*P23.^2 - P1.*P2.*P3) + ((P2.*P13 - P12.*P23).*(mu3 - ...
    z))./(P3.*P12.^2 - 2.*P12.*P13.*P23 + P2.*P13.^2 + P1.*P23.^2 - P1.* ...
    P2.*P3))));

xmin    =          -R;
xmax    =           R;
ymin    = @(x)     -sqrt(R^2 - x.^2);
ymax    = @(x)      sqrt(R^2 - x.^2);
zmin    = @(x,y)   -sqrt(R^2 - x.^2 - y.^2);
zmax    = @(x,y)    sqrt(R^2 - x.^2 - y.^2);

IPC  = integral3(funHndl,xmin,xmax,ymin,ymax,zmin,zmax, ...
                    "AbsTol",1e-10,"RelTol",1e-6);
end

