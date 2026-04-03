function smdLim = computeSmdLimGmm(gmmP,w,pp)
% propCovariance propagates covariance and yields limits for IPC and SMD
%
% INPUT: P       = [-] (3,3,N) Covariance matrix for each node.
%        pp      = [struct]    Postprocess structure.
%         
% OUTPUT: ipcLim = [-] (N,1) Limit IPC for each node
%         smdLim = [-] (N,1) Limit SMD for each node
%
% Documentation: https://docs.mosek.com/9.2/toolbox/tutorial-cqo-shared.html 
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
N      = pp.N;
n      = size(gmmP,4);
obj    = pp.obj;
ipcLim = pp.ipcLim;
HBR    = pp.HBR;
smdLim = nan(n,N);
for j = 1:n
    detP = nan(N,1);
    for i = 1:N
        detP(i)  = det(gmmP(:,:,i,j));                                     % [-] (N,1) Node-wise determinant of the positional covariance matrix
    end
    if strcmp(obj,'max_risk')    
        smdLim(j,:) = (sqrt(2)*HBR)^3./(3*exp(1)*sqrt(pi*detP)*ipcLim);         % [-] (N,1) SMD limit computed with the maximum risk formula
    elseif strcmp(obj,'risk')
        smdLim(j,:) = -2*log(3/w(j)*ipcLim*sqrt(pi*detP)/(sqrt(2)*HBR^3));           % [-] (N,1) SMD limit computed with Alfriend and Akella's formula applied to IPC
    else
        error('Invalid objective specification.');
    end
end
% If the IPClim is too high, the d_m^2 computed with these formulas can
% be negative. Since these nodes do not influence the maneuver, the SMD
% limit is set as very low in them.
smdLim(smdLim <= 0) = 1e-5;                                            % [-] (N,1) 

for j = 1:n
    IPC(j) = sqrt(2)*HBR^3/(3*sqrt(pi*detP(1)))*exp(-smdLim(j,1)'/2);
end
end