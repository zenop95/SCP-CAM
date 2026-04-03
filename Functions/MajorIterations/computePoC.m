function [PoC,smd,md,PoC_tot] = computePoC(relTraj,P,pp)
% computePoC computes the PoC and SMD for the tajectory
%
% INPUT: relTraj = [-] (3,N)   Relative trajectory 
%        P       = [-] (3,3,N) Positional ovariances
%        smdLim  = [-] (N,1)   Limits of the SMD per node
%        pp      = [struct]    Postprocess structure
%        
% OUTPUT: PoC      = [-]    (M,1) PoC for each encounter
%         smd      = [-]    (M,1) SMD for each encounter
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
M   = length(pp.secondary);
smd = nan(M,1);
PoC = nan(M,1);

for j = 1:M
    i   = pp.NCA(j);                                                       % [-] (1,1) Conjuntion node
    HBR = pp.secondary(j).HBR;                                             % [-] (1,1) Conjunction HBR
    r   = relTraj(1:3,i,j);                                                % [-] (3,1) Position part of the relative trajetory
    e2b2      = pp.e2b(:,:,j);     
    e2b2(2,:) = [];                                                        % [-] (2,3) B-plane transform that excludes the y component (eta of the B-plane)
    P_i       = e2b2*P(:,:,i,j)*e2b2';                                     % [-] (2,2) 2D covariance in the B-plane
    r         = e2b2*r;                                                    % [-] (2,1) Relative position in the B-plane 
    smd(j)    = dot(r,P_i\r);
    md(j)     = norm(r);
    switch lower(pp.PoCType)
        case 'maximum'
            PoC(j) = maximumPc(r,P_i,HBR)*pp.w(j);                         % [-] (1,1) Maximum PoC computed with Alfriend' and Akella's formula
        case 'constant'
            PoC(j) = constantPc(r,P_i,HBR)*pp.w(j);                        % [-] (1,1) PoC computed with Alfriend' and Akella's formula
        case 'chan'
            PoC(j) = poc_Chan(HBR,P_i,smd(j))*pp.w(j);                     % [-] (1,1) PoC computed with Chan's formula
        case 'alfano'
            PoC(j) = AlfanoPc(r,P_i,HBR)*pp.w(j);                          % [-] (1,1) PoC computed with Alfano's formula
        otherwise
            error('Invalid PoC type.')
    end
end
PoC(PoC>1) = 1;                                                            % [-] (N,1) To have nicer semilogy plots exclude undefined logarithms 
PoC_tot = PoCTot(PoC);
end