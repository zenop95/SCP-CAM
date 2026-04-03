function pp = runSCP(pp)
% runSCP
% Encapsulates the SCP pipeline: split GMM, set convex sets, find nodes,
% build AIDA assets, run major iterations, optional linear IPC rerun,
% optional SMD limit adaptation, validation, postprocessing.

% Optional preprocess for GMM order
if pp.gmmOrder > 1
    pp = splitInitial(pp);
end

% Define convexification sets / indicators
pp = defIndConv(pp);

% Firing nodes
pp = findFiringNodes(pp);

% ---- Major iterations
iterSmd = 0;
N       = pp.N;
newTraj = nan(6, N);

tic;
[newTraj, iterSmd, pp] = majorSolve(iterSmd, newTraj, pp);

% ---- Optional: Linear PoC constraint re-run
if pp.embedLinearConstraint
    pp.ipcConstr  = true;     % toggle linear constraint in majors
    prevMaxMin    = pp.iterMaxMin;
    pp.iterMaxMin = 1;        % single minor per major
    [~,~,pp]      = majorSolve(iterSmd, newTraj, pp);
    pp.iterMaxMin = prevMaxMin;
    if isfield(pp,'majorIter') && ~isempty(pp.majorIter)
        pp.t = pp.majorIter(end).t;
    end
end

% ---- Optional: SMD limit adaptation loop
if pp.adaptSmdLimit
    err       = 1;
    eta       = 0.1;             % stopping tolerance on relative error
    pp.nu_max = pp.nu_max / 10;   % tighten initial limit

    while err > eta && iterSmd < pp.iterMaxMaj
        if pp.fastEncounter
            pp = splitProbLimitShort(pp);
            [newTraj, iterSmd, pp] = majorSolve(iterSmd, newTraj, pp);
            prob = pp.majorIter(end).PcTot;
        else
            pp = splitProbLimitLong(pp);
            [newTraj, iterSmd, pp] = majorSolve(iterSmd, newTraj, pp);
            if isfield(pp.majorIter(end),'ipcTot')
                prob = max(pp.majorIter(end).ipcTot);
            else
                prob = pp.majorIter(end).PcTot;
            end
        end
        err = abs(prob - pp.lim) / max(pp.lim, eps);
    end
end

pp.simTime = toc;

end