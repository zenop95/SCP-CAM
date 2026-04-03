function pp = switchObj(pp)
%SwitchObj assigns the limit value to the relevant variable
% 
% INPUT: pp = [struct] postprocess structure
% 
% OUTPUT:  pp = [struct] postprocess structure
% 
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
switch pp.obj
    case 'miss_distance'
        pp.dMin     = pp.lim;
        pp.ipcLim   = nan;
    case 'max_risk'
        pp.ipcLim   = pp.lim;
    case 'risk'
        pp.ipcLim   = pp.lim;
end
end