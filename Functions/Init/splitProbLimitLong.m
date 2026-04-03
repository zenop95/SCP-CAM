function pp = splitProbLimitLong(pp)
% Adapt the IPC thresholds for the single conjunctions using a nonlinear
% program.
%
% INPUT: pp = [struct] Postprocess structure
%
% OUTPUT: pp = [struct] Postprocess structure
%
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

M          = length(pp.secondary);
n_c        = pp.n_conj;
lastIter   = pp.majorIter(end);
[ipcMax,i] = max(lastIter.ipcTot(2:end-1));
lim        = pp.lim;
if ipcMax < .99999*pp.lim
    cov     = lastIter.P;
    rN = nan(3,n_c); rO = rN; rHat = rN; rho0 = nan(n_c,1);
    for k = 1:n_c
        rN(:,k)   = lastIter.absP(:,i);
        rO(:,k)   = pp.ballisticTraj(1:3,i);
        rHat(:,k) = normalize(rN(:,k) - rO(:,k),'norm');
        rho0(k)   = norm(rN(:,k) - rO(:,k));
    end
%     alpha0  = -log10(pp.lims/pp.lim);
    alpha0  = -log10(normalize(lastIter.ipc(i,:)'/pp.lim,'norm',1));
    lb      = zeros(1,M+n_c);
    ub      = [3*ones(1,M) rho0'];
    alpha0(alpha0 > ub(1)) = ub(1);
    x0      = [alpha0; rho0];
    options = optimoptions( ...
                           'fmincon', ...
                           'Display',               'none',      ...
                           'ConstraintTolerance',   pp.lim/1e-6, ...
                           'Algorithm',             'sqp',       ...
                           'StepTolerance',         1e-10,       ...
                           'MaxIterations',         2e3,         ...
                           'MaxFunctionEvaluations',2e3          ...
                           ); 
    xnew    = fmincon(@(x)minFun(x),x0,[],[],[],[],lb,ub, ...
                @(x)nonLinConstr(x,pp.lim,rHat,rO,cov,pp),options);
    alpha = xnew(1:M);
    pp.lims = lim*10.^(-alpha);
end
end

function J = minFun(x)
    rho = x(end);
    J   = abs(rho);
end

function [c_in,c_eq] = nonLinConstr(x,lim,rHat,rO,cov,pp)
    n_c    = pp.n_conj;
    [~,i]  = max(pp.majorIter(end).ipcTot);
    order  = pp.gmmOrder;
    M      = n_c*order;
    lims   = nan(M,1); ipc = lims;
    alpha  = x(1:M);
    rho    = x(M+1:end);
    for k = 1:n_c
        for j = order*(k-1) + 1 : order*k
            HBR       = pp.secondary(j).HBR;
            r         = rho(k)*rHat(:,k) + rO(:,k) - pp.secondary(j).cart(1:3,i);                                                % [-] (3,1) Position part of the relative trajetory
            C         = cov(:,:,i,j);                                      % [-] (3,3) 3D covariance in the B-plane
            lims(j)   = lim*10^(-alpha(j));
            ipc(j)    = constantIpc(r,C,HBR)*pp.w(j);                      % [-] (1,1) PoC computed with Alfriend' and Akella's formula
        end
    end
    c_eq = lim - PoCTot(lims);
    c_in = ipc - lims;
end