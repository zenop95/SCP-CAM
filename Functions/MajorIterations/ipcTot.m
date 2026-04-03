function ipc_tot = ipcTot(ipc_single)
% PoCTot computes the total PoC when multiple conjunctions are considered.
% 
% INPUT: ipc_single = [-] (N,M) Vector with the PoC of single conjunctions
%        
% OUTPUT: ipc_tot   = [-] (N,1) PoC of all the conjunctions
% 
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

N = size(ipc_single,1);
M = size(ipc_single,2);
ipc_tot = nan(N,1);
for i = 1:N
    prob_none = 1;
    for j = 1:M
        prob_none = prob_none*(1-ipc_single(i,j));
    end
    ipc_tot(i) = 1 - prob_none;
end
end