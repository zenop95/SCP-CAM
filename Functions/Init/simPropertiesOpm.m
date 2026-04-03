function pp = simPropertiesOpm(opm,orbit,pp)
% simProperties initializes the simulation parameters. Modifications to the
% initialization need to be done inside the function.
%
% INPUT: orbit = string defining the orbit type as LEO or GEO
%
% OUTPUT: pp = structure with initialization parameters
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

%% Scaling factors and scaled initial variables
Lsc = 6378.173;
Vsc = sqrt(pp.mu/Lsc);
Tsc = Lsc/Vsc;
Asc = Vsc/Tsc;
pp.scaling     = [Lsc*ones(3,1); Vsc*ones(3,1); Asc*ones(3,1)];
pp.Tsc         = Tsc;

pp.HBR = pp.HBR/Lsc;
%% Fast encounter
pp.fastEncounter         = true;
pp.enableBplaneAvoidance = pp.fastEncounter*true;
pp.propCovariance        = false;

%% Time parameters
pp.cart_s = opm.state;
pp.x_s    = RV2MEE(pp.cart_s);
pp.x_d    = RV2MEE(pp.cart_s);
pp.uMax   = pp.opm.uMax;
pp.uMin   = 0;
pp.P0     = pp.cdm.cov1(1:3,1:3) + pp.cdm.cov2(1:3,1:3);
pp.dt     = pp.opm.dt;
pp.utc    = opm.epoch_spl; % utc at median simulation time t0
pp.order  = 2;
pp.N_forw = opm.N-1; % nodes after the median node
pp.N_back = 0; % nodes before the median node
pp.N      = pp.N_back + pp.N_forw+1;
pp.t0     = 0;                           %[s] median time where the covariance is known
pp.tb     = -pp.N_back*pp.dt;            %[s] starting time of simulation
pp.tf     = pp.t0+pp.N_forw*pp.dt;       %[s] ending time of simulation
pp.t      = linspace(pp.tb,pp.tf,pp.N);  %[s]
if pp.fastEncounter
    pp.NCA0   = pp.N_back+1; 
    pp.NCAf   = pp.N_back+1; % nodes in which CA constraint is active
else
    pp.NCA0   = 2; 
    pp.NCAf   = pp.N; % nodes in which CA constraint is active
end
pp.NSK0      = 2;
pp.NSKf      = pp.N; % nodes in which SK constraint is active
pp.et        = opm.et_start; %ephemeris time of simulation start
pp.dynamics  = 'aida'; %dynamics
pp.nDAVars   = 9;
pp.nDepVars  = 12;
%% Collision risk indicator selection
pp.lim       = 1e-6; %IPC limit
% pp.lim       = 2e3; %[m] miss distance limit
pp.obj       = 'risk';
pp.ipc_type  = 'constant'; if pp.fastEncounter; pp.ipc_type  = '2d'; end
pp           = switchObj(pp);


%% Optimization parameters  
pp.iterMaxMin      = 10; % Maximum number of minor iterations
pp.iterMaxMaj      = 10; % Maximum number of major iterations
pp.tolMin          = 1e-6; % [m] minor tolerance
pp.tolMaj          = 1e-3; % [m] major tolerance
pp.forwProp        = false; % Optimize using forward propagation or segmentation scheme
pp.enableTrAndVc   = true;
pp.trWeight        = 0;
pp.vcWeight        = 1e5;
pp.ctrlWeight      = pp.uMax;
pp.previousSolName = 'sol';

%% Toggle statation keping
pp.stationKeeping = false;  
pp.enableSkTarget = true;
pp.skSoft         = false;              % toggle station keeping as soft contraint
pp.skAlt          = false;
pp.skDev          = [0.05 0.05 1e4]';   % [deg] [deg] [m] station keeping lla box
pp.skTarget       = nan(6,1);           % Initialize variable
pp.loadTarget     = false;               % To avoid computing the target for the same scenario multiple times
pp.targWeight     = 100;                % Weight of the target SK objective
pp.skWeight       = 100;                % Weight of the slack latlon in objective for soft constraint
if strcmpi(orbit,'leo')
    pp.stationKeeping = false; 
    pp.skSoft         = false; 
end

%% Toggle smd maximum gradient constraint
pp.enableSmdGradConstraint = false;
pp.smdSoft                 = false;
if ~pp.enableSmdGradConstraint; pp.smdSoft = false; end
pp.smdGradWeight           = 0;  
pp.smdGradSlackWeight      = 1e-4;     % Weight of the slack grad in objective for soft constraint 
pp.maxIpcDeviation         = 0.1;   % [-]
pp.sdmGradTrig             = 0.1;  % [-] This parameter is fundamental because if chosen poorly it can lead to weird results with high DeltaV
pp.gradLimContinuation     = false;

%% Toggle minimum acceleration constraint
pp.enableHomotopy     = false; % still better not to use it
pp.homotopy = struct(...
                     'eps',        1e-6, ...
                     'betaTrig',   0.1 ...
                     );
pp.improveConvergence = 0;
pp.sl = pp.enableSkTarget*12;
pp.DvMax = pp.uMax*pp.scaling(7)*pp.dt*pp.Tsc;

end
