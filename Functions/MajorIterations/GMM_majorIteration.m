function [newTraj, majorIter, nu_max, breakFlag, pp] = ...
    majorIteration(t, majorIter, refute, oldTraj, nu_max, iter, k, G_db, pp)
% MAJORITERATION  One major-iteration cycle for trajectory optimisation.
%
%   Propagates primary/secondary trajectories, evaluates risk and covariance,
%   builds linearised maps, and calls the successive-convexification solver.
%
%   Inputs
%     t         (N,1)    Time vector
%     majorIter [struct] Previous major-iteration structure
%     refute    (1,1)    If true, reuse previously computed blocks
%     oldTraj   (6,N)    Relative trajectory from previous iteration
%     nu_max    (1,1)    Current trust-region parameter
%     iter      (1,1)    Major-iteration counter
%     k         scalar   Homotopy variable
%     G_db      scalar   Homotopy variable
%     pp        [struct] Parameter structure
%
%   Outputs
%     newTraj   (6,N)    Updated trajectory from linConvexSolve
%     majorIter [struct] Maps, costs, and diagnostics for this iteration
%     nu_max    (1,1)    Updated trust-region parameter
%     breakFlag (1,1)    True if solver reported infeasibility
%     pp        [struct] Updated parameter/postprocess structure
%
% Author: Zeno Pavanello
% E-mail: zpav176@aucklanduni.ac.nz
% Date:   2022-2024
%--------------------------------------------------------------------------


% %% ── 1. Build or reuse linearised maps ───────────────────────────────────
% if refute
%     [dynMaps, llaMaps, cartMaps, probMaps, state, uP, lla, constPos, ...
%      highRisk, gradOn, prob, P, smdLim, xi_dyn, xi_ca, gradLim] = ...
%         unpackMajorIter(majorIter);
% else
%     [dynMaps, llaMaps, cartMaps, probMaps, state, uP, lla, constPos, ...
%      highRisk, gradOn, prob, P, smdLim, xi_dyn, xi_ca, gradLim, C_p, pp] = ...
%         buildMaps(t, oldTraj, iter, pp);
% end
% 
% %% ── 2. Solve convex subproblem (minor iterations) ────────────────────────
% [newTraj, majorIter, breakFlag] = linConvexSolve( ...
%     dynMaps, llaMaps, cartMaps, probMaps, state, k, G_db, uP, lla, ...
%     constPos, prob, P, smdLim, highRisk, gradLim, gradOn, ...
%     xi_dyn, xi_ca, nu_max, pp);
% 
% if breakFlag
%     warning('majorIteration: infeasible solution — exiting early.');
%     return
% end
% 
% %% ── 3. Assemble major-iteration outputs ─────────────────────────────────
% majorIter = majBuild(majorIter, t, nu_max, prob, ...
%     dynMaps, llaMaps, cartMaps, probMaps, state, lla, constPos, C_p, P, ...
%     smdLim, highRisk, gradLim, gradOn, xi_dyn, xi_ca, pp);
% 
% end % majorIteration

function [newTraj, majorIter, nu_max, breakFlag, pp] = ...
    majorIteration(t, majorIter, refute, oldTraj, nu_max, iter, k, G_db, pp)

breakFlag = false;

if refute

    [dynMaps, llaMaps, cartMaps, probMaps, state, uP, lla, ...
     constPos, highRisk, gradOn, prob, P, smdLim, xi_dyn, ...
     xi_ca, gradLim] = unpackMajorIter(majorIter);

else

    [dynMaps, llaMaps, cartMaps, probMaps, state, uP, lla, ...
     constPos, highRisk, gradOn, prob, P, smdLim, xi_dyn, ...
     xi_ca, gradLim, C_p, pp] = buildMaps(t, oldTraj, iter, pp);

end

%% 2. Solve convex subproblem
% IMPORTANT:
% linConvexSolve must be updated to handle:
%   dynMaps : 9 x 6 x N x Np
%   state   : 6 x N x Np
%   P       : 6 x 6 x N x (Np*Ns)
% while retaining ONE shared primary control history.

[newTraj, majorIter, breakFlag] = linConvexSolve( ...
    dynMaps, llaMaps, cartMaps, probMaps, state, k, G_db, ...
    uP, lla, constPos, prob, P, smdLim, highRisk, gradLim, ...
    gradOn, xi_dyn, xi_ca, nu_max, pp);

if breakFlag
    warning('majorIteration:Infeasible', ...
            'Convex subproblem is infeasible.');
    return
end

%% 3. Assemble major-iteration outputs
majorIter = majBuild(majorIter, t, nu_max, prob, ...
    dynMaps, llaMaps, cartMaps, probMaps, state, lla, ...
    constPos, C_p, P, smdLim, highRisk, gradLim, ...
    gradOn, xi_dyn, xi_ca, pp);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [dynMaps, llaMaps, cartMaps, probMaps, state, uP, lla, constPos, ...
%           highRisk, gradOn, prob, P, smdLim, xi_dyn, xi_ca, gradLim, C_p, pp] = ...
%     buildMaps(t, oldTraj, iter, pp)
% % BUILDMAPS  Propagate trajectories and compute all linearised maps.
% M1 = length(pp.primary);
% M2 = length(pp.secondary);
% 
% %% 1a. Optional TCA refinement on first iteration (fast-encounter + GMM)
% % if iter == 1 && pp.fastEncounter && pp.toggleRefineTca && pp.gmmOrder > 1
% %     pp = refineTcaInit(oldTraj, iter, pp.N, M, t, pp);
% % end
% N = pp.N;
% 
% %% 1b. Propagate primary
% aidaInit(pp, 'primary', 0);
% [dynMaps, STM, llaMaps, cartMaps, state, cart, lla, xi_dyn, pp] = ...
%     propagateDA(iter, oldTraj, false, pp);
% 
% %pp.nomLon       = lla(2,1);
% %pp.primary.cart = state;
% 
% Np = numel(pp.primary);
% 
% for ip = 1:Np
%     pp.primary(ip).cart = cart(:,:,ip);
%     pp.primary(ip).state = state(:,:,ip);
% end
% 
% % Use a clearly defined reference primary for quantities that still
% % require a single nominal trajectory.
% lla = lla(:,:,1);
% pp.nomLon = lla(2,1);
% 
% % primaryX0 = ;
% % 
% % Np = size(primaryX0, 2);
% % 
% % cartP = nan(6, N, Np);
% % cartP(:,:,1) = cart;
% % 
% % for a = 2:Np
% %     dx0 = primaryX0(:,a) - primaryX0(:,1);
% % 
% %     for i = 1:N
% %         cartP(:,i,a) = cart(:,i) + STM(:,:,i)*dx0;
% %     end
% % end
% 
% %% 1c. Propagate secondaries (first iteration only)
% if iter == 1
%     [pp, NCA] = propagateSecondaries(state, M, N, pp);
%     pp.NCA    = NCA;
% end
% 
% %% 1d. Relative trajectories
% %relTraj = computeRelTraj(cart, M, N, pp);
% Ns = numel(pp.secondary);
% Np = numel(pp.primary);
% K  = Np * Ns;
% 
% relTraj = computeRelTraj(cart, Ns, N, pp);
% 
% if iter == 1
%     pp.lims = pp.lim * pp.w / pp.n_conj;
% end
% 
% uP = getControlHistory(iter, zeros(N,1), pp);
% 
% %% 1e. Covariance propagation
% [P, pp, C_p] = propCovariance(STM, cartMaps, N, true, pp);
% 
% %% 1f. B-plane transforms (fast encounter)
% pp = computeBplaneTransforms(cart, M, pp);
% 
% %% 1g. Risk metrics (PoC / IPC / SMD)
% r_exp = getReferenceRelTraj(iter, relTraj, pp);
% [prob, PoC, ipc, smd, highRisk] = computeRiskMetrics(relTraj, r_exp, P, pp);
% 
% %% 1h. SMD limit per node
% smdLim = computeSmdLim(P, relTraj(1:3,:,:), highRisk, pp);
% if strcmpi(pp.obj, 'miss_distance')
%     highRisk = smdLim > smd;
% end
% 
% %% 1i. Linear PoC maps
% [prob, probMaps, xi_ca] = computePocMaps(r_exp, P, prob, pp);
% 
% %% 1j. GEO station-keeping target (first iteration)
% if pp.enableSkTarget && iter == 1
%     pp = computeSkTarget(state, pp);
% end
% 
% %% 1k. SMD gradient-limit constraint
% [gradLim, gradOn, pp] = computeSmdGradConstraint(iter, smdLim, smd, ipc, N, pp);
% 
% %% 1l. Store ballistic baseline (first iteration)
% if iter == 1
%     pp = storeBallisticBaseline(cart, state, relTraj, ipc, smd, PoC, P, M, N, pp);
% end
% 
% %% 1m. Constant position term
% constPos = cart(1:3,:);
% 
% end % buildMaps

function [dynMaps, llaMaps, cartMaps, probMaps, state, uP, ...
          lla, constPos, highRisk, gradOn, prob, P, smdLim, ...
          xi_dyn, xi_ca, gradLim, C_p, pp] = ...
          buildMaps(t, oldTraj, iter, pp)

%% Dimensions
Np = numel(pp.primary);
Ns = numel(pp.secondary);
N  = pp.N;
K  = Np * Ns;

%% 1. Propagate all primary mixands
aidaInit(pp, 'primary', 0);

[dynMaps, STM, llaMaps, cartMaps, state, cart, ...
    lla_all, xi_dyn, pp] = propagateDA(iter, oldTraj, false, pp);

% Store each primary trajectory.
for ip = 1:Np
    pp.primary(ip).cart  = cart(:,:,ip);
    pp.primary(ip).state = state(:,:,ip);
end

% Select the nominal reference primary for legacy functions
% that still require one reference trajectory.
referenceIndex = 1;

pp.referencePrimaryIndex = referenceIndex;
pp.referencePrimaryCart  = cart(:,:,referenceIndex);
pp.referencePrimaryState = state(:,:,referenceIndex);

% Keep the complete geodetic propagation for all mixands.
pp.primaryLla = lla_all;

% Legacy single-reference LLA.
lla = lla_all(:,:,referenceIndex);
pp.nomLon = lla(2,1);

%% 2. Propagate all secondary mixands
if iter == 1

    referenceState = state(:,:,referenceIndex);

    [pp, NCA] = propagateSecondaries(referenceState, Ns, N, pp);

    pp.NCA = NCA;

end

%% 3. Construct primary-secondary pair information
% Pair ordering:
% k = (ip-1)*Ns + is

pairPrimary   = zeros(K,1);
pairSecondary = zeros(K,1);
pairWeight    = zeros(K,1);

relTraj = zeros(6,N,K);

for ip = 1:Np

    for is = 1:Ns

        pair = (ip-1)*Ns + is;

        pairPrimary(pair)   = ip;
        pairSecondary(pair) = is;

        % Component weights.
        wp = pp.primary(ip).w;
        ws = pp.secondary(is).w;

        pairWeight(pair) = wp * ws;

        % Relative Cartesian trajectory.
        relTraj(:,:,pair) = cart(:,:,ip) ...
                          - pp.secondary(is).cart;
    end
end

% Store pair metadata.
pp.pairPrimary   = pairPrimary;
pp.pairSecondary = pairSecondary;
pp.pairWeight    = pairWeight;

% Retain pair weights for risk aggregation.
pp.w = pairWeight;
pp.n_conj = K;

% Normalize defensively.
if sum(pp.w) <= 0
    error('buildMaps:InvalidWeights', ...
          'The primary-secondary pair weights sum to zero.');
end

pp.w = pp.w / sum(pp.w);

%% 4. Reference relative trajectories
if iter == 1
    r_exp = relTraj;
else
    r_exp = pp.majorIter(end).relP;

    if ~isequal(size(r_exp), size(relTraj))
        error('buildMaps:ReferenceTrajectorySize', ...
            ['Stored reference relative trajectories do not match ', ...
             'the current primary-secondary pair dimensions.']);
    end
end

%% 5. Common primary control history
uP = getControlHistory(iter, zeros(N,1), pp);

%% 6. Covariance propagation
% This function MUST be adapted to accept:
%   STM(:,:,i,ip)
% and the secondary covariance for each is.
%
% It must return pair-relative covariances:
%   P(:,:,i,pair)
%
% Do not pass the 4-D STM to the old scalar-primary implementation.
[P, pp, C_p] = propCovarianceGMM(STM, cartMaps, N, pp);

%% 7. B-plane transforms
% Long-term case: no B-plane transformation required.
if pp.fastEncounter
    pp = computeBplaneTransformsGMM(cart, Ns, pp);
else
    pp.e2b = [];
end

%% 8. Risk metrics
[prob, PoC, ipc, smd, highRisk] = ...
    computeRiskMetricsGMM(relTraj, r_exp, P, pp);

%% 9. SMD limits
smdLim = computeSmdLimGMM(P, relTraj(1:3,:,:), highRisk, pp);

if strcmpi(pp.obj, 'miss_distance')
    highRisk = smdLim > smd;
end

%% 10. Linear PoC maps
% Long-term IPC formulation does not use linear PoC constraints.
if pp.fastEncounter
    [prob, probMaps, xi_ca] = computePocMaps(r_exp, P, prob, pp);
else
    probMaps = [];
    xi_ca = [];
end

%% 11. Station-keeping target
if pp.enableSkTarget && iter == 1
    pp = computeSkTarget(state(:,:,referenceIndex), pp);
end

%% 12. SMD gradient constraint
[gradLim, gradOn, pp] = computeSmdGradConstraint( ...
    iter, smdLim, smd, ipc, N, pp);

%% 13. Store ballistic baseline
if iter == 1
    pp = storeBallisticBaselineGMM( ...
        cart, state, relTraj, ipc, smd, PoC, P, pp);
end

%% 14. Constant position term
% Preserve all primary mixands.
constPos = cart(1:3,:,:);

end





% ─────────────────────────────────────────────────────────────────────────

function pp = refineTcaInit(oldTraj, iter, N, M, t, pp)
% Propagate primary + secondaries before TCA refinement.
aidaInit(pp, 'primary', 0);
[~, ~, ~, ~, state, ~, ~, ~, pp] = propagateDA(iter, oldTraj, false, pp);
pp.primary.cart = state;

pp.w     = nan(M, 1);
pp.n_conj = M / pp.gmmOrder;
NCA      = [];
for j = 1:M
    aidaInit(pp, 'secondary', j);
    [pp.secondary(j), NCA] = propSecondary(pp.secondary(j), state, NCA, pp);
    pp.secondary(j).cart   = pp.secondary(j).x;
    for i = 1:N
        pp.secondary(j).x(:,i) = pp.cart2x(pp.secondary(j).x(:,i));
    end
    pp.w(j) = pp.secondary(j).w;
end
pp.NCA = NCA;
pp     = refineTca(M, t, state, pp);
end


% function [pp, NCA] = propagateSecondaries(state, M, N, pp) %#ok<INUSD>
% % Propagate all secondary objects on the first iteration.
% pp.w      = nan(M, 1);
% %pp.n_conj = M / pp.gmmOrder;
% pp.n_conj = K;
% NCA       = [];
% for j = 1:M
%     aidaInit(pp, 'secondary', j);
%     [pp.secondary(j), NCA] = propSecondary(pp.secondary(j), state, NCA, pp);
%     pp.secondary(j).cart   = pp.secondary(j).x;
%     pp.w(j)                = pp.secondary(j).w;
% end
% end

function [pp, NCA] = propagateSecondaries(referenceState, Ns, N, pp)

pp.secondaryWeights = zeros(Ns,1);
NCA = [];

for is = 1:Ns

    aidaInit(pp, 'secondary', is);

    [pp.secondary(is), NCA] = propSecondary( ...
        pp.secondary(is), referenceState, NCA, pp);

    pp.secondary(is).cart = pp.secondary(is).x;

    if ~isequal(size(pp.secondary(is).cart), [6,N])
        error('propagateSecondaries:TrajectorySize', ...
            'Secondary %d trajectory must be 6 x %d; received %s.', ...
            is, N, mat2str(size(pp.secondary(is).cart)));
    end

    pp.secondaryWeights(is) = pp.secondary(is).w;

end

end


% function relTraj = computeRelTraj(cart, M, N, pp)
% relTraj = nan(6, N, M);
% for j = 1:M
%     relTraj(:,:,j) = cart - pp.secondary(j).cart;
% end
% end


function relTraj = computeRelTraj(cart, Ns, N, pp)

Np = numel(pp.primary);
K  = Np * Ns;

relTraj = nan(6, N, K);

for ip = 1:Np
    for is = 1:Ns

        k = (ip-1)*Ns + is;

        relTraj(:,:,k) = cart(:,:,ip) - pp.secondary(is).cart;
    end
end
end

function uP = getControlHistory(iter, uP_default, pp)
if iter > 1
    uP = pp.majorIter(end).uP;
else
    uP = uP_default;
end
end


function pp = computeBplaneTransforms(cart, M, pp)
if pp.fastEncounter
    pp.e2b = nan(3, 3, M);
    for j = 1:M
        i = pp.NCA(j);
        pp.e2b(:,:,j) = eci2Bplane(cart(4:6,i), pp.secondary(j).cart(4:6,i));
    end
else
    pp.e2b = nan(1, 1, M);
end
end


function r_exp = getReferenceRelTraj(iter, relTraj, pp)
if iter == 1
    r_exp = relTraj;
else
    r_exp = pp.majorIter(end).relP;
end
end


function [prob, PoC, ipc, smd, highRisk] = computeRiskMetrics(relTraj, r_exp, P, pp) %#ok<INUSD>
if pp.fastEncounter
    [PoC, smd, prob] = computePoC(relTraj, P, pp);
    highRisk         = nan(0, size(relTraj,3));
    ipc              = nan;
else
    [ipc, smd, highRisk] = computeIpc(relTraj, P, pp);
    prob = ipcTot(ipc);
    PoC  = nan;
end
end


function [prob, probMaps, xi_ca] = computePocMaps(r_exp, P, prob, pp)
if pp.fastEncounter && ~strcmpi(pp.obj,'miss_distance') && pp.pocConstr
    [prob, probMaps, xi_ca] = pcDaMaps(r_exp, P, pp);
elseif ~strcmpi(pp.obj,'miss_distance') && pp.pocConstr
    error('Linear PoC constraint is not supported for long-term encounters.');
else
    probMaps = [];
    xi_ca    = [];
end
end


function pp = computeSkTarget(state, pp)
if strcmpi(pp.orbit, 'geo') && ~pp.loadTarget
    dt         = 5000 / pp.Tsc;
    tf         = pp.skLength * pp.T;
    cart_f     = state(:,end);
    [x0New, pp] = findGeoTarget(cart_f, tf, dt, pp.nomLon, pp);
elseif strcmpi(pp.orbit, 'geo') && pp.loadTarget
    x0New = pp.cart2x(load('skTarget.mat').target);
else
    x0New = state(1:6, end);   % LEO: return to nominal orbit
end
pp.skTarget   = x0New;
pp.cartTarget = pp.x2cart(x0New);
end


function [gradLim, gradOn, pp] = computeSmdGradConstraint(iter, smdLim, smd, ipc, N, pp)
gradLim      = 1e5 * ones(N, 1);
gradOn       = zeros(N, 1);
pp.gradScale = 1;
if pp.enableSmdGradConstraint && iter > 1
    [gradOn, gradLim, pp] = smdGradParams(iter, gradLim, smdLim, smd, ipc, pp);
end
pp.gradOn = gradOn;
end


function pp = storeBallisticBaseline(cart, state, relTraj, ipc, smd, PoC, P, M, N, pp)
pp.ballisticTraj    = cart;
pp.ballisticState   = state;
pp.ballisticIpc     = ipc;
pp.ballisticIpcTot  = ipcTot(ipc);
pp.ballisticSmd     = smd;
for j = 1:M
    pp.ballisticRelTraj(:,:,j) = cart - pp.secondary(j).cart;
end
for i = 1:N
    pp.ballisticSmdGrad(:,i) = smd_grad(P(:,:,i), relTraj(1:3,i,1));
end
if pp.fastEncounter
    pp.ballisticPoc    = PoC;
    pp.ballisticPoCTot = PoCTot(PoC);
end
end


function [dynMaps, llaMaps, cartMaps, probMaps, state, uP, lla, constPos, ...
          highRisk, gradOn, prob, P, smdLim, xi_dyn, xi_ca, gradLim] = ...
    unpackMajorIter(mi)
% UNPACKMAJORITERATION  Retrieve cached fields when refute == true.
dynMaps  = mi.dynMaps;   llaMaps  = mi.llaMaps;   cartMaps = mi.cartMaps;
probMaps = mi.probMaps;  state    = mi.state;      uP       = mi.uP;
lla      = mi.lla;       constPos = mi.constPos;   highRisk = mi.highRisk;
gradOn   = mi.gradOn;    prob     = mi.prob;        P        = mi.P;
smdLim   = mi.sqrMahaLim; xi_dyn = mi.xi_dyn;     xi_ca    = mi.xi_ca;
gradLim  = mi.gradLim;
end