function c = objective(m,pp)
% Objective computes the vector of linear coefficients to build the
% objective function of the optimization in MOSEK.
%
% INPUT: N           = [-] (1,1) Number of nodes in the problem. 
%        m           = [-] (1,1) Number of optimization variables per node
%        pp          = [struct]  Paramters structure.
% 
% OUTPUT: c          = [-] (mN,1) vector of the coefficients of the 
%                                 objective function 
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
ind  = pp.index;
N    = pp.N;
nOpt = m*N;
dt   = diff(pp.t); dt(end+1) = 0;
c    = zeros(m,N);                                                         % [-] (m,N) Initialize coefficients array
c(ind.ctrlCone,:) = dt/max(dt)*pp.ctrlWeight;                              % [-] (1,N) Cost coefficient of the control cone
c(ind.vcCone,:)   = pp.vcWeight;                                           % [-] (1,N) Cost coefficient of the virtual control cone
if pp.enableSmdGradConstraint && pp.smdSoft
        c(ind.gradSlack,:) = pp.smdGradSlackWeight;                        % [-] (3,N) Cost coefficient of the slack variables of the SMD gradient
end
if pp.stationKeeping && pp.skSoft
    c(ind.skSlack,:) = pp.skWeight;                                        % [-] (4,N) Cost coefficient of the slack variables of the SK
end
if pp.enableSkTarget
      targObj = pp.targWeight*ones(12,1);                                  % [-] (12,1) Cost coefficient of the slack variables of the target
else; targObj = []; 
end
c = [reshape(c,nOpt,1); targObj];                                          % [-] (nOpt,1) Reshape into vector
end