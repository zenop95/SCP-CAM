function dcm_BE = eci2Bplane(v_p,v_s)
% Eci2Bplane defines the B-plane for the short-term encounter
% INPUT: v_p = [-] (3,1) Velocity of the primary S/C in ECI coordinates
%        v_s = [-] (3,1) Velocity of the secondary S/C in ECI coordinates
%
% OUTPUT: dcm_BE = [-] (3,3) Direction cosine matrix from ECI to B-plane
%  
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

relVel   = v_p-v_s;
eta      = relVel/norm(relVel);
CrossVel = cross(v_s,v_p);
if norm(CrossVel) == 0     % particular case of miss distance = 0
   CrossVel = cross([1; 0; 0],eta);
   if norm(CrossVel) == 0     % particular case of miss distance = 0
    CrossVel = cross([0; 1; 0],eta);
   end
end
xi       = CrossVel/norm(CrossVel);
zeta     = cross(xi,eta);
zeta     = zeta/norm(zeta);
dcm_BE   = [xi'; eta'; zeta'];
end