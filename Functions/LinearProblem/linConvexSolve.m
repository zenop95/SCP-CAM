function [xNew, majorIter, flag] = linConvexSolve(dynMaps, llMaps, ...
    convMaps, probMaps, xOld, k, G_db, uP, lla, constPos, prob, P, ...
    smdLim, highRisk, gradLim, gradOn, xi_dyn, xi_ca, nu_max, pp)
% LINCONVEXSOLVE  Set up and solve the minor-iteration convex subproblem.
%
%   Runs successive minor iterations until convergence or the iteration
%   limit is reached, then returns the optimised state trajectory.
%
%   Inputs
%     dynMaps   (6,9,N)   Linear dynamics DA expansion
%     llMaps    (3,6,N)   Linear geodetic DA expansion
%     convMaps  (6,6,N)   Linear state-to-Cartesian DA expansion
%     xOld      (6,N)     Previous major-iteration solution
%     k         (1,1)     Homotopy parameter
%     G_db      (1,1)     Homotopy variable
%     uP        (N,1)     Previous homotopy variables
%     lla       (3,N)     Constant geodetic coordinates
%     constPos  (3,N)     Constant primary ECI position
%     prob      (1,1)     Collision probability
%     P         (3,3,N)   Relative-state covariance in ECI
%     smdLim    (N,1)     Node-wise SMD limit
%     highRisk  (N,1)     1 where SMD limit is violated, 0 otherwise
%     gradLim   (N,1)     Node-wise SMD gradient norm limit
%     gradOn    (N,1)     1 where SMD gradient constraint is active
%     xi_dyn    (6,N)     Non-linear parameter for trust region
%     xi_ca     (–)       Collision-avoidance trust-region parameter
%     nu_max    (1,1)     Maximum allowed NLI value
%     pp        struct    Parameter / postprocess structure
%
%   Outputs
%     xNew      (6,N)     Optimised state trajectory
%     majorIter struct    Contains all minor-iteration data
%     flag      (1,1)     True if solver returned an infeasible solution
%
% Author: Zeno Pavanello, 2022  |  zpav176@aucklanduni.ac.nz
% 2026: Code Improved using Claude (Sonnet 4.6)

%% ── 1. Initialise ────────────────────────────────────────────────────────
N       = pp.N;
M       = length(pp.secondary);
ind     = pp.index;
flag    = false;
majorIter = struct();

[x_exp, r_exp] = buildExpansionPoint(xOld, constPos, M, N, pp);
initEllP       = buildInitialEllipsoid(x_exp, P, smdLim, highRisk, M, N, pp);

%% ── 2. Minor-iteration loop ──────────────────────────────────────────────
iter    = 0;
err     = Inf;
tol     = pp.tolMin;
iterMax = pp.iterMaxMin;
skLim   = deg2rad(pp.skDev(1:2));

while err > tol && iter < iterMax
    iter = iter + 1;

    %% 2a. Solve convex subproblem
    [xNew, minorIter] = linConvexProblem(dynMaps, llMaps, convMaps, probMaps, ...
        x_exp, r_exp, xOld, constPos, prob, lla, initEllP, P, uP, smdLim, ...
        k, G_db, pp.skTarget, skLim, gradLim, gradOn, xi_dyn, xi_ca, nu_max, pp);

    %% 2b. Check feasibility — stop if infeasible
    if ~isFeasible(minorIter.feas); break; end

    %% 2c. Extract and post-process solution
    minorIter = extractSolution(xNew, minorIter, P, M, N, nu_max, ind, pp);
    xNew      = minorIter.x;

    %% 2d. Update ellipsoid expansion point for next iteration
    for j = 1:M
        initEllP(:,:,j) = minorIter.r_rel(:,:,j);
    end

    %% 2e. Compute convergence error
    err           = computeError(iter, minorIter, majorIter, pp);
    minorIter.err = err;

    %% 2f. Warn if virtual controls are non-negligible
    if sum(abs(minorIter.virtual)) > 1e-7
        warning('linConvexSolve: non-null virtual control = %.3e', ...
            sum(abs(minorIter.virtual)));
    end

    majorIter.minorIter(iter) = minorIter;
end

%% ── 3. Return result or flag infeasibility ───────────────────────────────
if isFeasible(minorIter.feas)
    xNew = minorIter.x;
else
    flag = true;
end

end % linConvexSolve


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_exp, r_exp] = buildExpansionPoint(xOld, constPos, M, N, pp)
% BUILDEXPANSIONPOINT  Choose expansion point from history or cold-start.
r_exp = nan(3, N, M);
if ~isfield(pp, 'majorIter')
    x_exp = [xOld; zeros(length(pp.index.ctrl), N)];
    for j = 1:M
        r_exp(:,:,j) = constPos - pp.secondary(j).cart(1:3,:);
    end
else
    prev = pp.majorIter(end);
    if pp.justInTime
        x_exp = [prev.x; prev.uNorm];   % (7,N) throttle as optimisation variable
    else
        x_exp = [prev.x; prev.u];       % (9,N) full control as optimisation variable
    end
    for j = 1:M
        r_exp(:,:,j) = prev.relP(:,:,j);
    end
end
end


function initEllP = buildInitialEllipsoid(x_exp, P, smdLim, highRisk, M, N, pp)
% BUILDINITIALELLIPSOID  Find the CA-constraint expansion point per node.
initEllP = nan(3, N, M);
for j = 1:M
    x_rel = nan(6, N);
    for i = 1:N
        x_rel(:,i) = pp.x2cart(x_exp(1:6,i)) - pp.secondary(j).cart(:,i);
    end
    r_rel = x_rel(1:3,:);

    if pp.fastEncounter
        N0 = pp.NCA(j);  Nf = pp.NCA(j);
    else
        N0 = pp.NCA0;    Nf = pp.NCAf;
    end

    [initEllP(:,:,j), ~] = initGuessEllipsoid(x_exp, r_rel, P(:,:,:,j), ...
        smdLim(:,j), highRisk(:,j), N0, Nf, pp.e2b(:,:,j), j, pp);
end
end


function minorIter = extractSolution(xNew, minorIter, P, M, N, nu_max, ind, pp)
% EXTRACTSOLUTION  Reshape solver output and compute derived quantities.
sl = pp.sl;
m  = pp.m;

if pp.enableSkTarget
    minorIter.slack = xNew(pp.index.targetSlack)';
end

xNew           = reshape(xNew(1:end-sl), m, N);
minorIter      = checkCones(xNew, minorIter, pp);
minorIter.x    = xNew(ind.state,:);

% Cartesian state
cart = nan(6, N);
for i = 1:N
    cart(:,i) = pp.x2cart(minorIter.x(:,i));
end
minorIter.p        = cart(1:3,:);
minorIter.v        = cart(4:6,:);
minorIter.ctrlCone = xNew(ind.ctrlCone,:);

% Relative position per secondary
for j = 1:M
    minorIter.r_rel(:,:,j) = minorIter.p - pp.secondary(j).cart(1:3,:);
end

% Control
if pp.justInTime
    minorIter.u = pp.uFix .* repmat(minorIter.ctrlCone, 3, 1);
else
    minorIter.u = xNew(ind.ctrl,:);
end

% Non-linearity indicator and homotopy variable
minorIter.xi = abs(minorIter.x - xNew(1:6,:)) * nu_max;
if pp.enableHomotopy
    minorIter.uP = xNew(m,:);
else
    minorIter.uP = nan;
end

% SMD gradient per node
for i = 1:N
    minorIter.smdGrad(:,i) = smd_grad(P(:,:,i), minorIter.r_rel(:,i));
end

minorIter.uNorm = normOfVec(minorIter.u);
end


function err = computeError(iter, minorIter, majorIter, pp)
% COMPUTEERROR  Convergence error: max change in relative position.
if iter > 1
    err = max(normOfVec(minorIter.r_rel - majorIter.minorIter(iter-1).r_rel));
elseif iter == 1 && isfield(pp, 'majorIter')
    err = max(normOfVec(minorIter.r_rel - pp.majorIter(end).minorIter(end).r_rel));
else
    err = Inf;
end
end


function ok = isFeasible(feasStr)
% ISFEASIBLE  True for accepted MOSEK feasibility statuses.
ok = strcmpi(feasStr, 'PRIMAL_AND_DUAL_FEASIBLE') || ...
     strcmpi(feasStr, 'UNKNOWN');
end