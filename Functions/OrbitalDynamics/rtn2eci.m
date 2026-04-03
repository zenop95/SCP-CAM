function [dcm,w] = rtn2eci(r, v)
% Calculates the cosine direction matrix to rotate the axis from the
% initial orbital reference frame LVLH to ECI.
% Usage:
%   dcm = eci2lvlh(r, v)
% Inputs:
%   - r: the orbital radius in ECI coordinates
%   - v: the spacecraft velocity in ECI coordinates
%   Note that r and v may be expressed in any measurement unit, as they
%   will be converted to versors; only their orientation is relevant.
% Outputs:
%   - dcm: a 3x3 matrix rotating LVLH into ECI; it might be employed as
%   follows: 
%           v_i = dcm_io * v_o
%   with v_i being a vector in the ECI frame, v_o in the orbital LVLH one.

%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
[dcm_o2e,w_o2e] = lvlh2eci(r,v);
dcm             = dcm_o2e*lvlh2rtn';
w               = lvlh2rtn*w_o2e;

end