function [ipc,smd,highRisk] = computeIpcGmm(gmmRelTraj,gmmP,w,smdLim,pp)
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
n           = size(gmmP,4);
N           = pp.N;
gmmSmd      = nan(n,N);
gmmHighRisk = zeros(n,N);
gmmIpc      = nan(n,N);
highRisk    = zeros(n,N);
for j = 1:n
    for i = 1:N
        r = gmmRelTraj(1:3,i,j);                                           % [-] (3,N)   Position part of the relative trajetory
        P_i = gmmP(:,:,i,j);                                                    % [-] (3,3,N) Positional combined covariance
        smd(j,i) = dot(r,P_i\r);
        switch pp.obj
            case 'max_risk'
                ipc(j,i) = maximumIpc(r,P_i,pp.HBR);                             % [-] (N,1) IPC computed with the maximum risk formula
            case 'risk'
                switch lower(pp.ipc_type)
                    case '2d'
                        ipc(j,i) = w(j)*pp.HBR^2/(2*sqrt(det(P_i)))*exp(-smd(j,i)/2);      % [-] (1,1) PC computed with Alfriend' and Akella's formula %%%%% da fixare nel multiple
                    case 'constant'
                        ipc(j,i) = w(j)*constantIpc(r,P_i,pp.HBR);                    % [-] (N,1) IPC computed with Alfriend' and Akella's formula 
                    case 'cuboid'
                        ipc(j,i) = w(j)*cuboidIpc(r,P_i,pp.HBR);                      % [-] (N,1) IPC computed with Zhang's cuboid formula
                    case 'numeric'
                        ipc(j,i) = w(j)*numericIpc(r,P_i,pp.HBR);                     % [-] (N,1) IPC computed with purely numerical integration
                end
        end
    smd(j,i) = w(j)*smd(j,i);
    end
end
highRisk(smd < smdLim) = 1;
ipc(ipc==0) = 1e-30;                                                       % [-] (N,1) To have nicer semilogy plots exclude undefined logarithms 
end