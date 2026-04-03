function pointEci = guessEllipseBplane(rEci,smdLim,PEci,e2b)
% pOnEllipsoidLine draws a segment from a point in the interior of the 
% ellipsoid to the ellipsoid border and finds the point of intersection 
% between this line and the ellipsoid border.
%
% INPUT: rEci   = [m]   relative distance between the two bodies expressed in 
%                       the ECI frame centered on the center of the ellipsoid.                    
%        smdLim = [m^2] objective value of the squared Mahalanobis dist.
%        PEci   = [m^2] covariance in ECI
%        pp     = []    postprocess structure
% OUTPUT: pointEci  = [m] 3x1 vector of the position of the point on the 
%                         ellipsoid border expressed in CW coordinates.
%         ell       = [ ] Ellipsoid structure.
%
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
%pass in 2d to solve B-plane problem
e2b(2,:)         = [];
rB               = e2b*rEci;
PB               = e2b*PEci*e2b';   
[semiaxes,cov2b] = defineEllipsoid(PB,smdLim);

%Define tot ellipse points
t          = 0:0.001:2*pi;
x          = semiaxes(1)*cos(t);
y          = semiaxes(2)*sin(t);
ellCov     = [x; y];
ellB       = nan(2,length(t));
for k = 1:length(t)
    ellB(:,k)   = cov2b*ellCov(:,k);
end

% find closest point to original
[~,b]    = min(normOfVec(ellB-rB));
pointB   = ellB(:,b);

% pass in 3D again
pointEci = e2b'*pointB;
end