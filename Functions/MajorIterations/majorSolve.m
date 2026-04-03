function [newTraj,iterf,pp] = majorSolve(iter0,newTraj,pp)
% MajorSolve performs the major-iteration logic that repeatedly calls the
%
% majorIteration solver, optionally updates trust-region (xi) parameters,
%
% evaluates convergence, and stores iteration data into pp.majorIter.
%
% INPUT:
%
%  iter0    = [-]      (1,1) Starting iteration number
%
%  newTraj  = [-]      (6,N) Initialized trajectory (state sequence)
%
%  pp       = [struct] Postprocess / options structure (uses fields such as
%
%                       N, iterMaxMaj, tolMaj, adaptTrustRegion, nus,
%
%                       nu_m, nu_M, nu_max, enableHomotopy, ctrlWeight, t, etc.)
%
% OUTPUT:
%
%  newTraj  = [-]      (6,N) Updated trajectory after major iterations
%
%  iterf    = [-]      (1,1) Final iteration number
%
%  pp       = [struct] Updated postprocess structure with pp.majorIter entries
%
% BEHAVIOR:
%
%  - Repeatedly calls majorIteration until convergence or limits are reached.
%
%  - If pp.adaptTrustRegion is true, calls updateXiMax to adapt nu_max and rho.
%
%  - Records per-iteration costs, control updates, and error measures in pp.
%
%  - May call homotopyUpdate when pp.enableHomotopy is true.
%
%  - Throws an error if majorIteration reports a breakFlag.
%
%  - Uses/updates fields pp.timeSubtr, pp.rho, pp.nus, and pp.majorIter.
%
% USAGE:
%
%  [newTraj, iterf, pp] = majorSolve(iter0, newTraj, pp);
%
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

N             = pp.N;
err           = 1; 
k             = 0;
nUpd          = 0;
G_db          = 1;
nu_max        = pp.nu_max;
t             = pp.t;
refute        = false;
majorIter     = struct();
it            = 0;
iter          = iter0;
while iter < pp.iterMaxMaj && it < pp.iterMaxMaj*2 && err > pp.tolMaj
    iter                 = iter + 1; 
    it                   = it + 1;
    if pp.adaptTrustRegion && it > 2 && ((pp.nus(it-1) == pp.nu_m && ...
       pp.nus(it-2) == pp.nu_m) || ...
       (pp.nus(it-1) == pp.nu_M && pp.nus(it-2) == pp.nu_M))
        warning('Minimum or Maximum \nu_{max} reached'); 
        break 
    end
    oldTraj              = reshape(newTraj,6*pp.N,1);
    [newTraj,majorIter,nu_max,breakFlag,pp] = ...
                             majorIteration(t,majorIter,refute, ...
                                    oldTraj,nu_max,iter,k,G_db,pp);
    if breakFlag
            error(['The optimizer yielded an infeasible formulation or there is ' ...
               'no need to optimize. Please check for any bugs present in ' ...
               'the propagation or in the convex problem formulation.']); 
    end
    if pp.adaptTrustRegion
        [nu_max,rho,linCost,nonLinCost,refute] = updateXiMax(pp,majorIter,nu_max,newTraj,iter);
        pp.timeSubtr = pp.timeSubtr + load("timeOut.dat");
        pp.rho(it) = rho;
        pp.nus(it) = nu_max;    
        if refute; iter = iter - 1; continue; end
    else
        nonLinCost = sum(normOfVec(majorIter.u))*pp.ctrlWeight;
        linCost = nonLinCost;
        refute = false;
    end
    pp.majorIter(iter)   = majorIter;
    pp.majorIter(iter).linCost = linCost;
    pp.majorIter(iter).J = nonLinCost;
    if iter > 1
        err = max(normOfVec(majorIter.u - pp.majorIter(iter-1).u));
%         err = abs(pp.majorIter(iter-1).J - linCost);
    end
    if iter0 > 0 && iter == iter0 + 1
        err = 1;
    end
    if pp.enableHomotopy
        [k,nUpd,G_db] = homotopyUpdate(pp,majorIter,k,nUpd,iter,hTrue);
    end
    pp.majorIter(iter).err = err;
    t = pp.majorIter(end).t;
end
iterf = iter;
end