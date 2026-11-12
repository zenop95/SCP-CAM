function pp = scenarioParams(up, dp)
% scenarioParams
% Builds the unscaled scenario: objects, nodes, time grid, toggles.
%
% Author: Zeno Pavanello
% E-mail: zpav176@aucklanduni.ac.nz
% Date:   2022-2024
%--------------------------------------------------------------------------

% Start from user params
pp = up;

% House defaults (for fields that user did not/sh8ould not set)
pp.mu               = dp.mu;
pp.flagRtn          = dp.flagRtn;
pp.dynamics         = dp.dynamics;
pp.nDAVars          = dp.nDAVars;
pp.nDepVars         = dp.nDepVars;
pp.adaptTrustRegion = dp.adaptTrustRegion;
pp.toggleRefineTca  = false; 

pp = generateInitGmmPS(pp); % GMM Splitting for Primary & Secondary objects

% Basic converters
pp.cart2x = @null;
pp.x2cart = @null;

% Discretization
nxOrb       = up.nxOrb;
pp.N_forw   = round(up.n_orbits_forward  * nxOrb);
pp.N_back   = round(up.n_orbits_backward * nxOrb);
pp.N        = pp.N_forw + pp.N_back + 1;

% Time grid (UNSCALED here, scaling later)
pp.dt_raw = pp.T / nxOrb;
pp.t0     = 0;
pp.tb     = -pp.N_back * pp.dt_raw;
pp.tf     =  pp.N_forw * pp.dt_raw;
pp.t      = linspace(pp.tb, pp.tf, pp.N)';

% Collision avoidance active window
pp.NCA0 = 2;
pp.NCAf = pp.N;

% Epoch handling (et computed after scaling)
pp.utc             = up.utc;
pp.timeSubtr       = 0;

% Objective and risk config
pp.obj      = up.obj;
pp.ipc_type = up.ipc_type;

% Linear PoC/IPC toggles
pp.pocConstr             = logical(up.pocConstr);

% Adapt SMD limit: auto-enable if there are multiple primaries and/or
% secondaries or GMM>1 for primary/secondary 
pp.adaptSmdLimit = logical((length(pp.primary) > 1 || length(pp.secondary) > 1 || pp.gmmOrder(1) > 1 || pp.gmmOrder(2) > 1) && ~pp.pocConstr);

% Optimizer settings
pp.iterMaxMaj = up.iterMaxMaj;
pp.iterMaxMin = up.iterMaxMin;
pp.tolMaj     = up.tolMaj;
pp.tolMin     = up.tolMin;

% Weights (ctrlWeight will be normalized as 1/N later)
pp.trWeight   = dp.trWeight;
pp.vcWeight   = dp.vcWeight;
pp.ctrlWeight = dp.ctrlWeight;

% Unscaled propulsion limits (scaled later)
pp.uMax_unscaled = up.uMax;
pp.uMin_unscaled = up.uMin;

% NLI limits (already scaled values; re-applied after scaling constants)
pp.nu_max = dp.nu_max;
pp.nu_m   = dp.nu_m;
pp.nu_M   = dp.nu_M;

% SMD gradient constraint (soft flag depends on enable flag)
pp.enableSmdGradConstraint = logical(up.enableSmdGradConstraint);
pp.smdSoft                 = pp.enableSmdGradConstraint && dp.smdSoft;
pp.smdGradSlackWeight      = dp.smdGradSlackWeight;
pp.maxIpcDeviation         = dp.maxIpcDeviation;
pp.sdmGradTrig             = dp.sdmGradTrig;
pp.gradLimContinuation     = dp.gradLimContinuation;

% Operational constraints
pp.justInTime    = dp.justInTime;
pp.uFixRTN       = dp.uFixRTN(:)/norm(dp.uFixRTN); % unit RTN direction
pp.canAlwaysFire = dp.canAlwaysFire;

% Homotopy
pp.enableHomotopy     = dp.enableHomotopy;
pp.homotopy           = dp.homotopy;
pp.improveConvergence = dp.improveConvergence;

% Bookkeeping
pp.nUpd = 0;

end