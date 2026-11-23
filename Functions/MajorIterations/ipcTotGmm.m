function ipc_totGmm = ipcTotGmm(ipcPair)
% ipcTotGmm computes the total IPoC when multiple conjunctions and mixands are considered.
% 
% INPUT: ipcPair = [-] (N,K) Vector with the IPoC of single conjunctions
% per mixand pair
%        
% OUTPUT: ipc_totGmm = [-] () IPoC of all the conjunctions
% 
% Author: 
% E-mail: 
%--------------------------------------------------------------------------
[N, K] = size(ipcPair); % K = Ns*Np
ipc_totGmm = nan(N,1);
for i = 1:N
    prob_none = 1;
    for j = 1:K
        prob_none = prob_none*(1-ipcPair(i,j));
    end
    ipc_totGmm(i) = 1 - prob_none;
end
end 