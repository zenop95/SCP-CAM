function pointEci = guessEllipseBplaneInPlane(x,rEci,smdLim,PEci,e2b)
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
%% ECI to RTN transformation
r2e  = rtn2eci(x(1:3),x(4:6));
% rRtn = r2e'*rEci;
% PRtn = r2e'*PEci*e2r';

%pass in 2d to solve B-plane problem
e2b(2,:)         = [];
r2b              = e2b*r2e;
nB               = r2b(:,3); % line made by the cut of the RT plane in the Bplane
nB               = nB/norm(nB);
rtB              = [0 1; -1 0]*nB;
rB               = e2b*rEci;
PB               = e2b*PEci*e2b';   
[semiaxes,cov2b] = defineEllipsoid(PB,smdLim);
rtCov            = cov2b'*rtB;
rtCovSc          = rtCov./semiaxes;
normrtCovSc      = norm(rtCovSc);
if normrtCovSc > 0
    pointCov = rtCovSc/normrtCovSc;
    pointB1  = cov2b*(pointCov.*semiaxes);
    pointB2  = cov2b*(-pointCov.*semiaxes);
    diff1 = norm(pointB1 - rB);
    diff2 = norm(pointB2 - rB);
    if diff1 >= diff2
        pointB = pointB2;
    else
        pointB = pointB1;
    end
else
    pointB    = cov2b*[semiaxes(1); 0];  %particular case of miss distance = 0: find point on semiminor axis
end
% pass in 3D again
pointEci = e2b'*pointB;
end