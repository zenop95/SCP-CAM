function poc_tot = PoCTot(poc_single)
% PoCTot computes the total PoC when multiple conjunctions are considered.
% 
% INPUT: poc_single = [-]      (M,1) Vector with the PoC of single conjunctions
%        s          = [struct] (M,1) Secondary structure
%        
% OUTPUT: poc_tot   = [-] (1,1) PoC of all the conjunctions
% 
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------


prob_none = 1;
M = length(poc_single);

for j = 1:M
    prob_none = prob_none*(1-poc_single(j));
end
poc_tot = 1 - prob_none;
end