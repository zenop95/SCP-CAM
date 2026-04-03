function [ipc,smd,highRisk] = computeIpc(relTraj,P,pp)
% computeIpc computes the IPC and SMD for the tajectory
%
% INPUT: relTraj = [-] (3,N)   Relative trajectory 
%        P       = [-] (3,3,N) Positional ovariances
%        smdLim  = [-] (N,1)   Limits of the SMD per node
%        pp      = [struct]    Postprocess structure
%        
% OUTPUT: ipc      = [-]    (N,1) IPC node-wise
%         smd      = [-]    (N,1) SMD node-wise
%         highRisk = [bool] (N,1) 1 where the IPC limit is violated                     
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
N        = pp.N;
M        = length(pp.secondary);
smd      = nan(N,M);
highRisk = zeros(N,M);
ipc      = zeros(N,M);

for j = 1:M
    HBR   = pp.secondary(j).HBR;
    for i = 1:N
        r   = relTraj(1:3,i,j);                                            % [-] (3,N)   Position part of the relative trajetory
        P_i = P(:,:,i,j);                                                  % [-] (3,3,N) Positional combined covariance
        smd(i,j) = dot(r,P_i\r);
            switch lower(pp.ipc_type)
                case 'maximum'
                    ipc(i,j) = maximumIpc(r,P_i,HBR)*pp.w(j);              % [-] (N,1) IPC computed with the maximum risk formula
                case 'constant'
                    ipc(i,j) = constantIpc(r,P_i,HBR)*pp.w(j);              % [-] (N,1) IPC computed with purely numerical integration
               case 'cuboid'
                    ipc(i,j) = cuboidIpc(r,P_i,HBR)*pp.w(j);               % [-] (N,1) IPC computed with Zhang's cuboid formula
                case 'serra'
                    ipc(i,j) = serraIpc(r,P_i,HBR,1e3)*pp.w(j);               % [-] (N,1) IPC computed with Zhang's cuboid formula
                case 'numeric'
                    ipc(i,j) = numericIpc(r,P_i,HBR)*pp.w(j);              % [-] (N,1) IPC computed with purely numerical integration
                otherwise
                    error('Invalid IPC type.')
            end
    end
    highRisk = ipc > pp.lim;
end
ipc(ipc>1) = 1;                                                            % [-] (N,1) Probability cannot be higher than 1 
end 