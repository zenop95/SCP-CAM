function [xNew, minorIter, scvx] = linConvexProblemAvA(p,x_exp,r_rel,P,smdLim,k,G_db, ...
    skTarget,skLim,gradLim,gradOn,xi_max,pp)
% linConvexProblem Solves the linear convex optimization problem of the
% multinode trajectory optimization using the conic optimization tool of 
% MOSEK.
%
% INPUT: dynMaps    = [-]    (6,9,N) Linear part of the dynamics DA expansion.
%        llMaps     = [-]    (3,6,N) Linear part of the geodetic DA expansion.
%        convMaps   = [-]    (6,6,N) Linear part of the conversion from state
%                                    to cartesian DA expansion.
%        x_exp      = [-]    (9,N)   Expansion point for linearizations.
%        x_const    = [-]    (9,N)   Constant part of the propagation.
%        pos_const  = [-]    (3,N)   Constant part of the primary ECI position.
%        ll_const   = [-]    (3,N)   Constant part of the geodetic coordinates.
%        r_rel      = [-]    (3,N)   Expansion point for CA constraint.
%        P          = [-]    (3,3,N) Covariance matrix of the relative state
%                                    expressed in ECI coordinates.
%        uP         = [-]    (N,1)   Previous iteration homotopy variables.
%        smdLim     = [-]    (N,1)   Node-wise limit of the SMD.
%        k          = [-]    (1,1)   Homotopy parameter.
%        highRisk   = [bool] (N,1)   0 if the SMD is respected, 1 otherwise.
%        skTarget   = [-]    (6,1)   Target state for final SK
%        skLim      = [-]    (2,1)   SK limits during the maneuver.
%        gradLim    = [-]    (N,1)   Limit value for the norm of SMD gradient
%        gradOn     = [bool] (N,1)   1 if SMD gradient is on, 0 otherwise.
%        nu         = [-]    (6,N)   Non-linear parameter for trust region.
%        xi_max     = [-]    (1,1)   Maximum allowed value for the NLI
%        pp         = [-]            Postprocess structure.
% 
% OUTPUT: xNew      = [-] (mN,1) Vector of the optimized variables.
%         minorIter = [struct]   Minor iteration structure.
%         scvx      = [struct]   SCVX structure with A matrices, limits,
%                                objective and cones
%
% Documentation: https://docs.mosek.com/9.2/toolbox/tutorial-cqo-shared.html 
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
%% Initialization of the problem
N       = pp.N;
m       = pp.m;
sl      = pp.sl;
nOpt    = m*N + sl;                                                        % [-] (1,1) Total length of the optimization vector
uMin    = pp.uMin/pp.uMax;                                                 % [-] (1,1) Minimum control
uMax    = 1;                                                               % [-] (1,1) Maximum control

NCA0    = pp.NCA0;                                                         % [-] (1,1) Starting node for CA
NCAf    = pp.NCAf;                                                         % [-] (1,1) Ending node for CA
scvx.A_rel   = []; scvx.A_ca     = []; scvx.b_rel    = [];
scvx.A_grad  = []; scvx.A_hom    = []; scvx.A_sk     = []; 
scvx.A_skTar = []; scvx.A_dx     = []; scvx.b_dx_up  = []; 
scvx.b_dx_lo = []; scvx.b_grad   = []; scvx.b_hom_lo = []; scvx.b_hom_up = [];
scvx.b_skTar = []; scvx.b_sk_lo  = []; scvx.b_sk_up  = []; 
scvx.A_alt   = []; scvx.b_alt_lo = []; scvx.b_alt_up = []; 
scvx.A_gradBound = []; scvx.b_gradBound = []; scvx.b_ca_lo = [];
scvx.limsUp = [];  scvx.limsLo = []; minorIter = [];
%% Limits and cost function
limsUp = zeros(nOpt,2); limsLo = zeros(nOpt,2); c = zeros(nOpt,2);
for j = 1:2
    ind = pp.index(j);
    [limsUp(:,j), limsLo(:,j)] = limits(x_exp(1:6,:,j),x_exp(1:6,:,2),N,m,...
                              uMax,gradLim,gradOn,p(j).nu,ind,xi_max,pp);  % [struct] Define the limits for the optimization variables
    c(:,j)      = objective(N,m,ind,j,pp);                                   % [struct] Define the coefficients for the objective function
end    
scvx.limsUp = sum(limsUp,2);
scvx.limsLo = sum(limsLo,2);
scvx.c      = sum(c,2);
%% Basic problem
[A_dx1,b_dx_up1,b_dx_lo1] = ...
                 nliConstr(N,m,nOpt,x_exp(:,:,1),p(1).nu,pp.index(1));     % [struct] Build the NLI constraints
[A_dx2,b_dx_up2,b_dx_lo2] = ...
                 nliConstr(N,m,nOpt,x_exp(:,:,2),p(2).nu,pp.index(2));     
scvx.A_dx    = [A_dx1; A_dx2];
scvx.b_dx_up = [b_dx_up1; b_dx_up2];
scvx.b_dx_lo = [b_dx_lo1; b_dx_lo2];

% if pp.justInTime
%     [scvx.A_dyn, scvx.b_dyn] = contConstrJit(N,m,sl,dynMaps,x_const, ...
%                                                              x_exp,pp);  % [struct] Build the dynamics constraints in the Just-in-time scenario
% else
x_e1      = reshape([x_exp(:,:,1); zeros(m-9,N)],m*N,1); zeros(sl,1);      % [-] (mN,1) Reshape expansion point for propagation
x_e2      = reshape([zeros(m/2,N); x_exp(:,:,2); ...
                            zeros(m/2-9,N)],m*N,1); zeros(sl,1);           % [-] (mN,1) Reshape expansion point for propagation
[A_dyn1,b_dyn1] = contConstr(N,m,sl,p(1).dynMaps,p(1).state, ...
                                x_e1,pp.index(1),pp);                      % [struct] Build the dynamics constraints
[A_dyn2,b_dyn2] = contConstr(N,m,sl,p(2).dynMaps,p(2).state, ...
                                x_e2,pp.index(2),pp);                      % [struct] Build the dynamics constraints
scvx.A_dyn = [A_dyn1; A_dyn2];
scvx.b_dyn = [b_dyn1; b_dyn2];
% end

[scvx.A_rel,scvx.b_rel] = relPosConstrAvA(N,m,nOpt,p,x_exp,pp);                              % [struct] Build the relative position constraints
if pp.fastEncounter; N0 = pp.NCA; Nf = pp.NCA; 
else; N0 = pp.NCA0; Nf = pp.NCAf; end
[A,b,minorIter] = collAvoidConstr(m,nOpt,N0,Nf, ...
   r_rel,P,smdLim,pp.e2b,1,minorIter,pp);                                  % [struct] Build the CA constraints
scvx.A_ca    = [scvx.A_ca; A];
scvx.b_ca_lo = [scvx.b_ca_lo; b];
%% Additional constraint
% if pp.enableSmdGradConstraint
%     [scvx.A_grad, scvx.b_grad] = ...
%                    smdGradConstr(m,nOpt,NCA0,NCAf,gradOn,P,scvx.A_grad,pp);% [struct] Build the SMD grad constraints
%     [scvx.A_gradBound, scvx.b_gradBound] = smdGradBound(nOpt,NCA0, ...
%                                                  NCAf,gradOn,m,gradLim,pp);
% end
% 
% if pp.stationKeeping
%    [scvx.A_sk,scvx.b_sk_up,scvx.b_sk_lo] = ...
%      skConstr(m,nOpt,skLim,ll_const(1:2,:),llaMaps(1:2,:,:),x_exp,pp);     % [struct] Build the SK constraints
% end
% 
% if pp.altSk
%    [scvx.A_alt,scvx.b_alt_up,scvx.b_alt_lo] = ...
%      altConstr(m,nOpt,ll_const(3,:),llaMaps(3,:,:),x_exp,pp);       % [struct] Build the SK constraints
% end
% 
% if pp.enableSkTarget
%    [scvx.A_skTar,scvx.b_skTar] = targConstr(m,N,nOpt,skTarget,pp);         % [struct] Build the targeting constraints
% end
% 
% if pp.enableHomotopy
%     [scvx.A_hom,scvx.b_hom_lo,scvx.b_hom_up] = ...
%                        homConstr(N,m,sl,uP,k,G_db,uMin,uMax,pp);           % [struct] Build the homotopy constraints
% end

%% Create MOSEK problem structure
prob = createProblem(m,N,scvx,pp);

%% Run optimization
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-11;
param.MSK_DPAR_INTPNT_TOL_PFEAS    = 1e-11;
[~,res]         = mosekopt('minimize echo(0)',prob,param);
% param.MSK_IPAR_INFEAS_REPORT_AUTO = 1;
% [~,res]         = mosekopt('minimize anapro',prob,param);
try
    minorIter.feas  = res.sol.itr.prosta;                                  % [str] Extract state of the optimization solution
    if strcmpi(minorIter.feas,'unknown')
        warning('The optimizer could not find an optimal solution')
    end
    xNew            = res.sol.itr.xx;                                      % [-] (mN,1) Extract optimized vector
%     mosekResponse(res);
catch; error(['error in MOSEK: ', res.rcodestr, '. ', res.rmsg]); 
end
end