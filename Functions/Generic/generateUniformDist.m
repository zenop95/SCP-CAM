function set = generateUniformDist(lims,n)
%GENERATEUNIFORMDIST Summary of this function goes here
%   Detailed explanation goes here
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

set = unifrnd(lims(1),lims(2),[1 n]);
end

