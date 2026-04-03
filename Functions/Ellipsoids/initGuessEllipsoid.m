function [new_r_rel, ellipsoids] = ...
                  initGuessEllipsoid(x_exp,r_rel,P,smdLim,highRisk,NCA0,NCAf,e2b,j,pp)
% pOnEllipsoid draws a segment from a point in the interior of the 
% ellipsoid to the ellipsoid border and finds the point of intersection 
% between this line and the ellipsoid border.
%
% INPUT: r_rel    = [-]    (3xN) Distance of the primary from the secondary
%                   expressed in ECI coordinates.                    
%        P        = [-]    (3x3xN) Covariance matrices of the relative 
%                   distance expressed in ECI coordinates.
%        smdLim   = [-]    (Nx1) Limit SMD value.
%        highRisk = [bool] (Nx1) 0 if the SMD is respected, 1 otherwise.
%        NCA0     = [-]    (1x1) Starting node for CA constraint.
%        NCAf     = [-]    (1x1) Ending node for CA constraint.
%        pp       = [struct] Postprocess structure.
%
% OUTPUT: new_r_rel  = [-] (3xN) array with the shifted position of the
%                      original points inside the ellipsoid
%         ellipsoids = [struct] Ellipsoid structure.
%
% Validated using the script testExpansionPoint
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
% Initialize vector of the guessed points on the ellipsoids surfaces
ellipsoids = struct('semiaxes', nan,      ...
                    'Dr',       nan(3,1), ...
                    'z',        nan(3,1), ...
                    'cov2CW',   nan(3,3)  ...
                    );
%% Find new points outside the ellipsoid
new_r_rel = r_rel;
if pp.enableBplaneAvoidance
    i = pp.NCA(j);
    new_r_rel(:,i) = ...
        guessEllipseBplane(r_rel(:,i),smdLim,P(:,:,i),e2b);                % Find the point on the ellipse for short-term encounters          
%         guessEllipseBplaneInPlane(x_exp(:,i),r_rel(:,i),smdLim,P(:,:,i),e2b);                % Find the point on the ellipse for short-term encounters          
else
    for i = NCA0:NCAf 
        if highRisk(i)
            [new_r_rel(:,i),ellipsoids] = ...
                guessEllipsoid(r_rel(:,i),P(:,:,i),smdLim(i));             % Find the point on the ellipsoid for long-term encounters    
        end
    end
end
end