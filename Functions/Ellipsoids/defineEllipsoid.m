function [semiaxes,cov2eci] = defineEllipsoid(P,smdLim)
% DefineEllipsoid Finds the semiaxes of the ellipsoid and the transfomration 
% between the original frame and the covariance frame. It also works in 2D
% for the B-plane ellipse.
%
% INPUT: P      = [-] (3x3) Covariance matrix of the relative distance
%                     expresse in ECI coordinates
%        smdLim = [-] (1x1) Limit value of the square Mahalanobis dist
% 
% OUTPUT: semiaxes  = [-] (3x1) Semiaxes of the keep-out ellipsoid.
%         cov2eci   = [-] (3x3) DCM rotating from covariance frame to ECI
% 
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
% try chol(P);
% catch; error('The covariance matrix is not positive definite'); end
% Compute the rotation from the CW ref system to a covariance ref system
[V,D]   = eig(P);
[a,b]   = sort(diag(D),'descend');
V       = V(:,b);
D       = diag(a);
cov2eci = V;
if det(cov2eci) < 0; cov2eci(:,1) = -cov2eci(:,1); end

% Normalized semiaxes of the ellipsoid corresponding to smdLim
a        = sqrt(smdLim*D(1,1));
b        = sqrt(smdLim*D(2,2));
semiaxes = [a; b];                                                         % In B-plane ellipse case only two semiaxes
if length(D) == 3
    semiaxes = [semiaxes; sqrt(smdLim*D(3,3))];                            % In ellipsoid case add the third semiaxis
end
end