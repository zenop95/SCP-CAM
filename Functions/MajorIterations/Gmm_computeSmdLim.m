function smdLim = Gmm_computeSmdLim(P,pp) 
% INPUT:  P = [-] (3,3,N,K)
% OUTPUT: smdLim = [-] (N,K) Limit SMD per mixand pair for each node
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
N      = pp.N;
K     = size(P,4);

obj    = pp.obj;
ipcLim = pp.ipcLim;
HBR    = pp.HBR;

smdLim = nan(N, K);

for j = 1:K
    detP = nan(N,1);
    for i = 1:N
        % Pi = P(:,:,i,j);
        % Pi = (Pi+Pi')/2;
        detP(i)  = det(P(:,:,i,j));                                     % [-] (N,1) Node-wise determinant of the positional covariance matrix
    end
    %if strcmp(obj,'max_risk')    
    %   smdLim(j,:) = (sqrt(2)*HBR)^3./(3*exp(1)*sqrt(pi*detP)*ipcLim);         % [-] (N,1) SMD limit computed with the maximum risk formula
    if strcmp(obj,'risk')
        smdLim(:,j) = -2*log(3/w(j)*ipcLim*sqrt(pi*detP)/(sqrt(2)*HBR^3));           % [-] (N,1) SMD limit computed with Alfriend and Akella's formula applied to IPC
    else
        error('Invalid objective specification.');
    end
end
% If the IPClim is too high, the d_m^2 computed with these formulas can
% be negative. Since these nodes do not influence the maneuver, the SMD
% limit is set as very low in them. [A negative threshold means the risk inequality
% is already satisfied for any nonnegative SMD.]
smdLim(smdLim <= 0) = 1e-5;                                            % [-] (N,1) 

% for j = 1:K
%     IPC(j) = sqrt(2)*HBR^3/(3*sqrt(pi*detP(1)))*exp(-smdLim(1,j)'/2);
% end
end
