function up = userParams()
% USERPARAMS Defines specific parameters for the algorithm. The code should
% be carefully edited by the user.
% 
% Author: Zeno Pavanello
% E-mail: zpav176@aucklanduni.ac.nz
% Date:   2022-2024
%--------------------------------------------------------------------------

%% Scenario selection
up.orbit        = 'leo';           % 'leo' | 'geo'
up.scenario     = 'single-short';   % 'single-short' | 'single-long' | 'multiple-short' | 'repeating-short'
up.ESACase      = 1;               % '1-2170', only used if up.scenario == 'single-short'
up.gmmOrder     = 1;              % Must be odd number

%% discretization per orbit
up.nxOrb = 60;

%% Time window [in number of orbits from first TCA]
up.n_orbits_forward  = 0;
up.n_orbits_backward = 1;

%% Propulsion (UNSCALED)
up.uMax = 3.75e-7;   % [km/s^2]
up.uMin = 0;      % [km/s^2]

%% Constraints & toggles
up.pocConstr               = false;  % if true => minor iters forced to 1
up.enableSmdGradConstraint = false;

%% Optimizer limits
up.iterMaxMaj = 10;
up.iterMaxMin = 10;
up.tolMaj     = 1e-3;
up.tolMin     = 1e-6;

%% objective/risk settings
up.obj      = 'miss_distance';          % 'risk' | 'max_risk' | 'miss_distance' | '2d'
up.PoCType  = 'Chan';          % 'Constant' | 'Maximum' | 'Chan' | 'Alfano'
up.ipc_type = 'Constant';      % 'Constant' | 'Serra' | 'Maximum' | 'Cuboid'
up.mdLim    = 2;               % [km] 
up.pocLim   = 1e-6;

%% Station-keeping
up.stationKeeping = false;
up.skSoft         = false;
up.altSk          = false;   
up.enableSkTarget = true;

%% Epoch at TCA (UTC string)
up.utc = '06/03/2020 15:00:00';

%% AIDA environment models
up.aida.flag1    = 0;    % atmosphere flag (1:non-rotating, 2:rotating)
up.aida.flag2    = 0;    % SRP flag (1:no shadow, 2:Earth cylindrical shadow, 3:Earth biconical shadow, 4:Earth and Moon cylindrical shadow, 5:Earth biconical and Moon cylindrical shadow, 6:Earth and Moon biconical shadow)
up.aida.flag3    = 0;    % third body flag (1:Moon, 2:Moon and Sun)
up.aida.gravOrd  = 0;    % Order of the EGM2008 gravitational model

end