function [k,nUpd,G_db] = homotopyUpdate(pp,majorIter,k,nUpd,iter,hTrue)
% homotopyUpdate updates the homotopy parameter k based on 
% the formula presented in ref.
%
% INPUT: pp    = parameters structure
%        iter  = major iteration number
%
% OUTPUT: k    = [] homotopy parameter
%         nUpd = [] number of times the homotopy parameter has been updated
%
% Reference: Malyuta, D., & Acikmese, B. (2021). Fast Homotopy for 
% Spacecraft Rendezvous Trajectory Optimization with Discrete Logic. 
% http://arxiv.org/abs/2107.07001
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

% Update parameters
delta0    = 1e-1;
delta1    = 1e-2;
eps       = 1e-2;
betaTrig  = pp.homotopy.betaTrig;
rho       = delta1/delta0;
nMax      = pp.iterMaxMaj;
% Update rule
% if iter > 1
%     del = abs(sum(normOfVec(pp.majorIter(end).u))-sum(normOfVec(majorIter.u)));
% elseif iter == 1 && isfield(pp,"majorIter")
%     del = abs(sum(normOfVec(pp.majorIter(end).u))-sum(normOfVec(majorIter.u)));
% else; del = 0;
% end
% if  del <= betaTrig && any(hTrue)
    k    = log(1/eps-1)/(rho^(nUpd/nMax)*delta0);
    nUpd = nUpd + 1;
% end
u_db = pp.homotopy.u_db;
gMax = 1 - pp.uMin/pp.uMax;
G_db = -k*exp(k*(u_db)/gMax)/(gMax*(exp(k*(u_db)/gMax)+1)^2)*u_db + 1 - ...
        1/(1+exp(k*(u_db)/gMax)) + 1/(1+exp(k*gMax));
end