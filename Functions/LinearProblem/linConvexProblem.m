function [xNew, minorIter, scvx] = linConvexProblem(dynMaps, llaMaps, ...
    convMaps, probMaps, x_exp, r_exp, x_const, pos_const, probab, ll_const, ...
    r_rel, P, uP, smdLim, k, G_db, skTarget, skLim, gradLim, gradOn, ...
    xi_dyn, xi_ca, nu_max, pp)
% LINCONVEXPROBLEM  Build and solve the linearised convex subproblem via MOSEK.
%
%   Assembles all constraint matrices (dynamics, collision avoidance,
%   station-keeping, homotopy, etc.), creates the MOSEK problem structure,
%   and returns the optimised variable vector.
%
%   Inputs
%     dynMaps   (6,9,N)   Linear dynamics DA expansion
%     llaMaps   (3,6,N)   Linear geodetic DA expansion
%     convMaps  (6,6,N)   Linear state-to-Cartesian DA expansion
%     x_exp     (9,N)     Linearisation expansion point
%     r_exp     (3,N)     CA-constraint expansion point
%     x_const   (9,N)     Constant part of propagation
%     pos_const (3,N)     Constant primary ECI position
%     probab    (1,1)     Collision probability
%     ll_const  (3,N)     Constant geodetic coordinates
%     r_rel     (3,N,M)   Relative-position expansion point
%     P         (3,3,N,M) Relative-state covariance in ECI
%     uP        (N,1)     Previous homotopy variables
%     smdLim    (N,1)     Node-wise SMD limit
%     k         (1,1)     Homotopy parameter
%     G_db      (1,1)     Homotopy variable
%     skTarget  (6,1)     Final station-keeping target state
%     skLim     (2,1)     Station-keeping limits during manoeuvre
%     gradLim   (N,1)     Node-wise SMD gradient norm limit
%     gradOn    (N,1)     1 where SMD gradient constraint is active
%     xi_dyn    (6,N)     Non-linear trust-region parameter
%     xi_ca     (–)       CA trust-region parameter
%     nu_max    (1,1)     Maximum allowed NLI value
%     pp        struct    Parameter / postprocess structure
%
%   Outputs
%     xNew      (mN,1)  Optimised variable vector
%     minorIter struct  Minor-iteration diagnostics
%     scvx      struct  A-matrices, limits, objective, and cone definitions
%
% Author: Zeno Pavanello, 2022  |  zpav176@aucklanduni.ac.nz
% 2026: Code Improved using Claude (Sonnet 4.6)

%% ── 1. Initialise ────────────────────────────────────────────────────────
N    = pp.N;
m    = pp.m;
sl   = pp.sl;
M    = length(pp.secondary);
nOpt = m*N + sl;

uMin = pp.uMin / pp.uMax;
uMax = 1;

scvx      = initScvxStruct();
minorIter = struct();

%% ── 2. Limits and objective ──────────────────────────────────────────────
[scvx.limsUp, scvx.limsLo] = limits(x_exp(1:6,:), r_exp, N, m, ...
    uMax, gradLim, gradOn, xi_dyn, xi_ca, nu_max, pp);
scvx.c = objective(m, pp);

%% ── 3. Dynamics continuity constraint ───────────────────────────────────
if pp.justInTime
    [scvx.A_dyn, scvx.b_dyn] = contConstrJit(N, m, sl, dynMaps, x_const, x_exp, pp);
else
    [scvx.A_dyn, scvx.b_dyn] = contConstr(N, m, sl, dynMaps, x_const, x_exp, pp);
end
minorIter.dynResidual = scvx.b_dyn;

%% ── 4. Collision avoidance constraint ───────────────────────────────────
[scvx.A_ca, scvx.b_ca_lo, scvx.b_ca_up, minorIter] = ...
    buildCaConstraint(M, m, nOpt, r_rel, P, smdLim, probMaps, pos_const, ...
                      r_exp, probab, x_exp, minorIter, pp);

%% ── 5. Optional constraints ──────────────────────────────────────────────
if pp.enableSmdGradConstraint
    [scvx.A_grad, scvx.b_grad] = smdGradConstr(m, nOpt, pp.NCA0, pp.NCAf, ...
        gradOn, P, scvx.A_grad, pp);
    [scvx.A_gradBound, scvx.b_gradBound] = smdGradBound(nOpt, pp.NCA0, ...
        pp.NCAf, gradOn, m, gradLim, pp);
end

if pp.stationKeeping
    [scvx.A_sk, scvx.b_sk_up, scvx.b_sk_lo] = skConstr(m, nOpt, skLim, ...
        ll_const(1:2,:), llaMaps(1:2,:,:), x_exp, pp);
end

if pp.altSk
    [scvx.A_alt, scvx.b_alt_up, scvx.b_alt_lo] = altConstr(m, nOpt, ...
        ll_const(3,:), llaMaps(3,:,:), x_exp, pp);
end

if pp.enableSkTarget
    [scvx.A_skTar, scvx.b_skTar] = targConstr(m, N, nOpt, skTarget, pp);
end

if pp.enableHomotopy
    [scvx.A_hom, scvx.b_hom_lo, scvx.b_hom_up] = homConstr(N, m, sl, ...
        uP, k, G_db, uMin, uMax, pp);
end

%% ── 6. Solve with MOSEK ──────────────────────────────────────────────────
prob = createProblem(m, N, scvx, pp);

param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-11;
param.MSK_DPAR_INTPNT_TOL_PFEAS    = 1e-11;
[~, res] = mosekopt('minimize echo(0)', prob, param);

[xNew, minorIter] = extractMosekResult(res, minorIter);

end % linConvexProblem


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scvx = initScvxStruct()
% INITSCVXSTRUCT  Return an empty scvx struct with all expected fields.
fields = {'A_rel','A_ca','b_rel','A_grad','A_hom','A_sk','A_skTar', ...
          'A_dx','b_dx_up','b_dx_lo','b_grad','b_hom_lo','b_hom_up', ...
          'b_skTar','b_sk_lo','b_sk_up','A_alt','b_alt_lo','b_alt_up', ...
          'A_gradBound','b_gradBound','b_ca_lo','b_ca_up'};
scvx = cell2struct(repmat({[]}, numel(fields), 1), fields);
end


function [A_ca, b_ca_lo, b_ca_up, minorIter] = buildCaConstraint( ...
    M, m, nOpt, r_rel, P, smdLim, probMaps, pos_const, r_exp, probab, ...
    x_exp, minorIter, pp)
% BUILDCACONSTRAINT  Build collision-avoidance constraint matrices.
A_ca = []; b_ca_lo = []; b_ca_up = [];

if pp.pocConstr && pp.fastEncounter
    if strcmpi(pp.obj, 'miss_distance')
        [A_ca, b_ca_up, b_ca_lo] = linMdConstr(M, m, nOpt, probMaps, pos_const, r_exp, pp);
    else
        [A_ca, b_ca_up, b_ca_lo] = linPoCConstr(M, m, nOpt, probMaps, probab, x_exp, pp);
    end
    return
end

for j = 1:M
    if pp.fastEncounter
        N0 = pp.NCA(j);  Nf = pp.NCA(j);
    else
        N0 = pp.NCA0;    Nf = pp.NCAf;
    end
    [A, b_up, b_lo, minorIter] = collAvoidConstr(m, nOpt, N0, Nf, ...
        r_rel(:,:,j), P(:,:,:,j), smdLim(:,j), pp.e2b(:,:,j), j, minorIter, pp);
    A_ca    = [A_ca;    A   ];
    b_ca_lo = [b_ca_lo; b_lo];
    b_ca_up = [b_ca_up; b_up];
end
end


function [xNew, minorIter] = extractMosekResult(res, minorIter)
% EXTRACTMOSEKRESULT  Parse MOSEK output and handle errors.
try
    minorIter.feas = res.sol.itr.prosta;
    if strcmpi(minorIter.feas, 'UNKNOWN')
        warning('linConvexProblem: MOSEK could not find an optimal solution.');
    end
    xNew = res.sol.itr.xx;
catch
    error('linConvexProblem: MOSEK error — %s. %s', res.rcodestr, res.rmsg);
end
end