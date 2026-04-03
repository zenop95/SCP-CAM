function pp = splitProbLimitShort(pp)
% Adapt the PoC thresholds for the single conjunctions using a nonlinear
% program.
%
% INPUT: pp = [struct] Postprocess structure
%
% OUTPUT: pp = [struct] Postprocess structure
%
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

M     = length(pp.secondary);
n_c   = pp.n_conj;
lastIter = pp.majorIter(end);
lim     = pp.lim;
cov     = lastIter.P;
rN = nan(3,n_c); rO = rN; rHat = rN; rho0 = nan(n_c,1); beta0 = zeros(n_c,1);
for k = 1:n_c
    start     = pp.gmmOrder*(k-1) + 1;
    finish    = start + pp.gmmOrder - 1;
    i         = pp.NCA(start);
    rN(:,k)   = lastIter.absP(:,i);
    rO(:,k)   = pp.ballisticTraj(1:3,i);
    rHat(:,k) = normalize(rN(:,k) - rO(:,k),'norm');
    rho0(k)   = norm(rN(:,k) - rO(:,k));
    limConj   = pp.lims(start:finish);
    probConj  = lastIter.pc(start:finish);
    for j = 1:pp.gmmOrder
        beta0(k)  = beta0(k) + (abs(limConj(j) - probConj(j))/limConj(j) < 1/10);
    end
end

alpha0  = -log10(normalize(lastIter.pc/pp.lim,'norm',1));
% alpha0  = -log10(normalize(pp.lims/pp.lim,'norm',1));
x0      = [alpha0; rho0];
lb      = zeros(1,M+n_c);
ub      = [3*ones(1,M) rho0'];
options = optimoptions( ...
                       'fmincon',                           ...
                       'Display',               'none',     ...
                       'ConstraintTolerance',   pp.lim/1e6, ...
                       'Algorithm',             'sqp',      ...
                       'StepTolerance',         1e-10,      ...
                       'MaxIterations',         2e3,        ...
                       'MaxFunctionEvaluations',2e3         ...
                       );
xnew    = fmincon(@(x)minFun(x,beta0),x0,[],[],[],[],lb,ub, ...
            @(x)nonLinConstr(x,pp.lim,rHat,rO,cov,pp),options);
alpha = xnew(1:M);
rho = xnew(M+1:end);
pp.lims = lim*10.^(-alpha);
end

function J = minFun(x,beta0)
    n_c  = length(beta0);
    rho  = x(end - n_c + 1: end);
    J    = dot(rho,beta0);
end

function [c_in,c_eq] = nonLinConstr(x,lim,rHat,rO,cov,pp)
    n_c   = pp.n_conj;
    order = pp.gmmOrder;
    M     = n_c*order;
    lims  = nan(M,1); 
    pc    = nan(M,1);
    alpha = x(1:M);
    rho   = x(M+1:end);
    for k = 1:n_c
        for j = order*(k-1) + 1 : order*k
            HBR       = pp.secondary(j).HBR;
            i         = pp.NCA(j);
            r         = rho(k)*rHat(:,k) + rO(:,k) - pp.secondary(j).cart(1:3,i);                                                % [-] (3,1) Position part of the relative trajetory
            e2b       = pp.e2b([1 3],:,j);     
            C         = e2b*cov(:,:,i,j)*e2b';                                 % [-] (2,2) 2D covariance in the B-plane
            rB        = e2b*r; 
            lims(j)   = lim*10^(-alpha(j));                                         % [-] (1,1) SMD limit computed with Maximum Alfriend and Akella's formula applied to PC
            pc(j)     = constantPc(rB,C,HBR)*pp.w(j);                          % [-] (1,1) PoC computed with Alfriend' and Akella's formula
        end
    end
    c_eq = lim - PoCTot(lims);
    c_in = pc - lims;
end