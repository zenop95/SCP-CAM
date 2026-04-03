function maps = dat2maps(datFile,nDAVars,nDepVars,N)
% Dat2maps reads the dat file created by DACE and turns it into linear maps
%
% INPUT: datFile  = string with the name of the .dat file
%        nDAVars  = number of independent DA variables
%        nDepVars = number of dependent variables
%        N        = number of nodes in the optimization
%
% OUTPUT: maps = [nDepVars x nDAVars x N] array of the linear maps
%
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
    DAtab = readtable(datFile, 'Delimiter', 'tab');
    DAmaps = table2array(DAtab);
    DAmaps(isnan(DAmaps(:,1)),:) = [];
    mapsTr = reshape(DAmaps',[nDAVars,nDepVars,N]); 
    maps = nan(nDepVars,nDAVars,N);
    for i = 1:N; maps(:,:,i) = mapsTr(:,:,i)'; end
end