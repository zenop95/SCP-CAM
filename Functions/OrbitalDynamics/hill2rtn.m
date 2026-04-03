function dcm = hill2rtn()
%HILL2RTN DCM that transforms the Hill frame into the RTN frame

% Zeno Pavanello 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

dcm = lvlh2rtn()*lvlh2hill()';
end
