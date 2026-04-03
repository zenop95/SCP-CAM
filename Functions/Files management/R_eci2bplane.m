function R = R_eci2bplane(v_1_eci,v_2_eci) 
% Compute the rotation matrix from ECI to Bombardelli-defined Bplane
% -------------------------------------------------------------------------
% INPUTS:
% - v_1_eci: primary object velocity in ECI reference frame
% - v_2_eci: secondary object velocity in ECI reference frame
% -------------------------------------------------------------------------
% OUTPUT:
% - R : rotation matrix from ECI to Bplane
% -------------------------------------------------------------------------
% Author:        Marco Felice Montaruli, Politecnico di Milano, 27 November 2020
%                e-mail: marcofelice.montaruli@polimi.it


    u_ksi = cross(v_2_eci,v_1_eci)/norm(cross(v_2_eci,v_1_eci));
    u_eta = (v_1_eci - v_2_eci)/norm(v_1_eci - v_2_eci);
    u_zeta = cross(u_ksi,u_eta);
    R = [u_ksi'; u_eta'; u_zeta'];

end