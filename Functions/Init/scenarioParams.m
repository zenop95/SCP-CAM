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

% House defaults (for fields that user did not/should not set)
pp.mu               = dp.mu;
pp.flagRtn          = dp.flagRtn;
pp.dynamics         = dp.dynamics;
pp.nDAVars          = dp.nDAVars;
pp.nDepVars         = dp.nDepVars;
pp.adaptTrustRegion = dp.adaptTrustRegion;
up.singleObject     = true;
pp.toggleRefineTca  = false; 

% Initial orbit & objects
if strcmpi(up.orbit,'geo')
    pp = generateInitLongGeo(pp);
else
    switch lower(pp.scenario)
        case 'single-short'
            pp = generateInitShort(pp);
        case 'single-long'
            pp = generateInitLong(pp);
        case 'multiple-short' 
            pp = generateInitMultipleStarlink(pp);
        case 'repeating-short'    
            pp = generateInitMultipleGmm(pp);
        otherwise
            error('pp.scenario must be set to a valid entry')
    end
end
% Encounter mode (affects windows later)
pp.enableBplaneAvoidance = pp.fastEncounter;

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
if pp.fastEncounter
    pp.NCA0 = pp.N_back + 1;
    pp.NCAf = pp.N_back + 1;
else
    pp.NCA0 = 2;
    pp.NCAf = pp.N;
end

% Station-keeping active window
pp.NSK0 = 2;
pp.NSKf = pp.N;

% Epoch handling (et computed after scaling)
pp.utc             = up.utc;
pp.timeSubtr       = 0;

% Objective and risk config
pp.obj      = up.obj;
pp.PoCType  = up.PoCType;
pp.ipc_type = up.ipc_type;

% Linear PoC/IPC toggles
pp.pocConstr             = logical(up.pocConstr);

% Adapt SMD limit: auto-enable if there are multiple secondaries or GMM>1
pp.adaptSmdLimit = logical((length(pp.secondary) > 1 || ...
                            pp.gmmOrder > 1) && ~pp.pocConstr);

% Station-keeping intent
pp.stationKeeping = logical(up.stationKeeping);
pp.skSoft         = logical(up.skSoft);

% SK target enabled only if there is a forward horizon
pp.enableSkTarget = up.enableSkTarget*(up.n_orbits_forward > 0);
pp.sl             = pp.enableSkTarget*12;
pp.skTarget       = nan(6,1);           
pp.loadTarget     = dp.loadTarget;           
pp.skLength       = dp.skLength;
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

% SK details (scale later)
pp.skDev      = dp.skDev;   % [deg]
pp.altLimKm   = dp.altLimKm;   % [km]
pp.targWeight = dp.targWeight;
pp.skWeight   = dp.skWeight;

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