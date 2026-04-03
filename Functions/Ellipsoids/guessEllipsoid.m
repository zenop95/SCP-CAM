function [pointEci, ell] = guessEllipsoid(Dr,P,smdLim)
% pOnEllipsoidLine draws a segment from a point in the interior of the 
% ellipsoid to the ellipsoid border and finds the point of intersection 
% between this line and the ellipsoid border.
%
% INPUT: Dr         = [m] relative distance between the two bodies
%                         expressed in the CW (or B-space)frame centered 
%                         on the center of the ellipsoid.                    
%        P          = [m^2] covariance matrix of the relative distance
%                           expressed in CW (or B-space) coordinates.
%        sqrMahaObj = [m^2] objective value of the squared Mahalanobis dist.
%
% OUTPUT: pointCW   = [m] 3x1 vector of the position of the point on the 
%                         ellipsoid border expressed in CW coordinates.
%         ell       = [ ] Ellipsoid structure.
%
% Validated using the script testExpansionPoint
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
[semiaxes,cov2eci] = defineEllipsoid(P,smdLim);
DrCov              = (cov2eci'*Dr)./semiaxes;
% Point of intersection between the surface and the line connecting the
% center of the ellipsoid and the point Dr
pointCov  = DrCov/norm(DrCov);
pointEci = cov2eci*(pointCov.*semiaxes);

ell = struct('semiaxes', semiaxes, ...
             'Dr',       Dr, ...
             'z',        pointEci, ...
             'cov2CW',   cov2eci ...
             );
end