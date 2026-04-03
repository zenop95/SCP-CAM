function dp = defaultParams()
% defaultParams
% House defaults that typical users should not edit directly.
%
% Author: Zeno Pavanello
% E-mail: zpav176@aucklanduni.ac.nz
% Date:   2022-2024
%--------------------------------------------------------------------------

dp.mu = 398600.4418;

% frames, flags
dp.flagRtn   = true;           % thrust in RTN if true

% weights (ctrlWeight will be set as 1/N later)
dp.trWeight   = 0;
dp.vcWeight   = 1e6;
dp.ctrlWeight = 1; % placeholder; actual used = 1/N

% caching/previous solutions
dp.loadPreviousSol = false;
dp.previousSolName = 'sol';
dp.adaptTrustRegion = false;

% dynamics backend
dp.dynamics  = 'aida';
dp.nDAVars   = 9;
dp.nDepVars  = 12;

% NLI limits (scaled)
dp.nu_max = 1e-3;
dp.nu_m   = 1e-10;
dp.nu_M   = 1e-2;

% station-keeping defaults
dp.skDev          = [0.05 0.05]';       % [deg] lat/lon bounds (soft/hard downstream)
dp.altLimKm       = [35783 35790.7]';   % [km] GEO box; scaled later
dp.targWeight     = 1e3;                % soft target weight
dp.skWeight       = 100;                % soft geodetic slack weight
dp.skLength       = 14;                 % Days of GEO SK validity after end of maneuver
dp.loadTarget     = false;              

% SMD gradient constraint defaults
dp.enableSmdGradConstraint = false;
dp.smdSoft            = false;   % enabled only when enableSmdGradConstraint=true
dp.smdGradSlackWeight = 1e-1;
dp.maxIpcDeviation    = 0.1;
dp.sdmGradTrig        = 0.1;
dp.gradLimContinuation = false;

% operational constraints
dp.justInTime    = false;
dp.uFixRTN       = [0;0;1]; % unit direction in RTN when JIT active
dp.canAlwaysFire = true;

% homotopy/min-accel convergence helper (disabled by default)
dp.enableHomotopy = false;
dp.homotopy = struct( ...
    'eps',      1e-2, ...
    'betaTrig', 0.1,  ...
    'u_db',     1e-5  ...
);

% misc
dp.improveConvergence = 0;

end