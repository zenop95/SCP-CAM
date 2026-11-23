function [ipcPair,smd,highRisk] = Gmm_computeIpc(relTraj,P,wPair,smdLim,pp)
% Gmm_computeIpc computes the IPC and SMD for the trajectory
%
% INPUT: relTraj = [-] (3,N,K)   Relative trajectory 
%        P       = [-] (3,3,N,K) Positional ovariances
%        wPair   = [-] (Kx1) 
%        smdLim  = [-] (N,K)   Limits of the SMD per node
%        pp      = [struct]    Postprocess structure
%        
% OUTPUT: ipcPair  = [-]    (N,K) IPC node-wise (already weighted)
%         smd      = [-]    (N,K) SMD node-wise
%         highRisk = [bool] (N,K) 1 where the IPC limit is violated                     
%         
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
N           = pp.N;
K           = size(relTraj,3);

wPair = wPair(:);
smd         = zeros(N,K);
ipcPair     = zeros(N,K);
highRisk    = zeros(N,K);

for j = 1:K
    for i = 1:N
        r = relTraj(1:3,i,j);                                           % [-] (3,N)   Position part of the relative trajetory
        P_i = P(:,:,i,j);
        Pi = (Pi+Pi')/2;
        [~,flag] = chol(Pi);
        if flag ~= 0
            error('Covariance is not positive definite.');
        end

       
        smd(i,j) = r'*(Pi\r);
        switch lower(pp.obj)
            % case 'max_risk'
            %    ipcPair(i,j) = maximumIpc(r,Pi,pp.HBR);                             % [-] (N,1) IPC computed with the maximum risk formula
            case 'risk'
                switch lower(pp.ipc_type)
                    % case '2d'
                    %    ipcPair(i,j) = wPair(j)*pp.HBR^2/(2*sqrt(det(Pi)))*exp(-smd(i,j)/2);      % [-] (1,1) PC computed with Alfriend' and Akella's formula %%%%% da fixare nel multiple
                    case 'constant'
                        ipcPair(i,j) = wPair(j)*constantIpc(r,Pi,pp.HBR);                    % [-] (N,1) IPC computed with Alfriend' and Akella's formula 
                    % case 'cuboid'
                    %    ipcPair(i,j) = wPair(j)*cuboidIpc(r,Pi,pp.HBR);                      % [-] (N,1) IPC computed with Zhang's cuboid formula
                    % case 'numeric'
                    %    ipcPair(i,j) = wPair(j)*numericIpc(r,Pi,pp.HBR);                     % [-] (N,1) IPC computed with purely numerical integration
                end
        end
    end
end
highRisk(smd < smdLim) = 1;
ipc(ipc==0) = 1e-30;                                                       % [-] (N,1) To have nicer semilogy plots exclude undefined logarithms 
end