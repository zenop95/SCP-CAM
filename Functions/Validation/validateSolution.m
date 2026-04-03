function [pp] = validateSolution(pp)
%Validate performs the validation propagation of the trajectory
% 
% INPUT:  pp = postprocess structure
%       
% OUTPUT: pp = postprocess structure
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
N                                  = pp.N;                                 % [-] (1,1) Number of nodes in the optimization
M                                  = length(pp.secondary);                 % [-] (1,1) Number of secondary objects
t                                  = pp.t;                   % [-] (1,1) Time step
pp.u                               = pp.majorIter(end).u;                  % [-] (3,N) Optimezed control history
pp.uRtn                            = pp.majorIter(end).uRtn;               % [-] (3,N) Optimezed control history
pp.uEci                            = pp.majorIter(end).uEci;               % [-] (3,N) Optimezed control history
pp.u(abs(pp.u) < pp.uMin/pp.uMax)  = 0;                                    % [-] (3,N) Exclude controls lower than minimum threshold
pp.uNorm                           = normOfVec(pp.u);                      % [-] (1,N) Node-wise norm of the control histroy
[~,STM,~,cartMaps,~,cart,lla,~,pp] = propagateDA(1,nan,true,pp);           % [-]       Validation propagation
pp.validationAbsTraj               = cart;                                 % [-] (6,N) Validated cartesian absolute trajectory
for j = 1:M
    pp.validationTraj(:,:,j)       = cart - pp.secondary(j).cart;          % [-] (6,N) Validated cartesian relative trajectory
end
pp.lla                             = lla;                                  % [-] (6,N) Validated cartesian trajectory

%% Propagate the covariance
P  = propCovariance(STM,cartMaps,pp.N,true,pp);                        % [-] (3,3,N) Propagated covariance
smdLim = pp.majorIter(end).sqrMahaLim;%computeSmdLim(P,pp.validationTraj(1:3,:),pp.majorIter(1).highRisk,pp);                                              % [-] (N,1) Node-wise SMD limit

%% Rotate DeltaV and acceleration to show it in LVLH frame
pp.uLvlh = zeros(3,N);                                                     % [-] (3,N) Acceleration in LVLH frame
r2eP     = nan(3,3,N);
% w = [0; -pp.n; 0];
for i = 1:N
    for j = 1:M
        l2eS          = lvlh2eci(pp.secondary(j).cart(1:3,i), ...
                                   pp.secondary(j).cart(4:6,i));           % [-] (3,3) DCM from secondary LVLH to ECI
        pp.pLvlh(:,i,j) = l2eS'*pp.validationTraj(1:3,i,j);                    % [-] (3,N,M) Relative position in secondary LVLH frame
    end
    if pp.uNorm(i) ~= 0
        r2eP(:,:,i)   = rtn2eci(pp.validationAbsTraj(1:3,i), ...
                                   pp.validationAbsTraj(4:6,i));           % [-] (3,3) DCM from primary LVLH to ECI
        pp.uEci(:,i) = r2eP(:,:,i)'*pp.u(:,i);                             % [-] (3,N) Acceleration in primary LVLH frame
    end
end
dt = diff(t); dt(end+1) = dt(end);
pp.dv     = pp.uRtn.*repmat(dt',3,1)*pp.uMax*pp.scaling(4)*1000;                        % [m/s] (3,N) Node-wise Delta V components in LVLH
pp.DvNorm = normOfVec(pp.dv);                                              % [m/s] (1,N) Node-wise Delta V magnitude
pp.DvTot  = sum(pp.DvNorm);                                                % [m/s] (1,1) Total Delta V

%% Propagate after validation to check target before and after maneuver
pp = validateGeo(pp);
%% Probability of Collision
if pp.enableBplaneAvoidance
    a = pp.PoCType;
    pp.PoCType   = 'Maximum';                                             
    [pp.PoCConst,pp.md,~,pp.PoCTotMaximum] = computePoC(pp.validationTraj,P,pp);
    pp.PoCType   = 'Constant';                                             
    [pp.PoCConst,pp.md,~,pp.PoCTotConstant] = computePoC(pp.validationTraj,P,pp);
    pp.PoCType   = 'Chan';                                             
    [pp.PoCChan,pp.md,~,pp.PoCTotChan]      = computePoC(pp.validationTraj,P,pp);
    pp.PoCType   = 'Alfano';                                             
    [pp.PoCAlfano,~,pp.md,pp.PoCTotAlfano]  = computePoC(pp.validationTraj,P,pp);
    pp.PoCType   = a;
    pp.md = pp.md*pp.scaling(1);
    smd          = nan;
    ipc                  = nan;
    ipc_tot      = nan;
else
    [ipc,smd] = computeIpc(pp.validationTraj,P,pp);                 % [-] (N,1) Node-wise IPC and SMD
    ipc_tot    = ipcTot(ipc);
    % pp.Pc = nan;
    % pp.PcTot = sum(ipc);                                                   % [-] (1,1) Probability of collision in the long-term encounter case
    % x1 = pp.ballisticTraj(:,1).*pp.scaling(1:6)*1e3;
    % x2 = pp.secondary.cart(:,1).*pp.scaling(1:6)*1e3;
    % D = diag(pp.scaling(1:6)*1e3);
    % C1 = D*pp.majorIter(1).covPrim(:,:,1)*D;
    % C2 = D*pp.secondary.covariance(:,:,1)*D;    
    % pp.PoCHallBall = Pc3D_Hall(x1(1:3),x1(4:6),C1,x2(1:3),x2(4:6),C2,pp.secondary.HBR*pp.scaling(1)*1e3);
    % x1 = pp.validationAbsTraj(:,2).*pp.scaling(1:6)*1e3;
    % x2 = pp.secondary.cart(:,2).*pp.scaling(1:6)*1e3;
    % D = diag(pp.scaling(1:6)*1e3);
    % C1 = D*pp.majorIter(end).covPrim(:,:,2)*D;
    % C2 = D*pp.secondary.covariance(:,:,2)*D;    
    % pp.PoCHallMan = Pc3D_Hall(x1(1:3),x1(4:6),C1,x2(1:3),x2(4:6),C2,pp.secondary.HBR*pp.scaling(1)*1e3);
end
% x1 = pp.validatedAbsTraj(:,1).*pp.scaling(1:6)*1e3;
% D = diag(pp.scaling(1:6));
% C1 = D*pp.majorIter(1).covPrim(:,:,1)*D;
% C2 = D*pp.secondary.covariance(:,:,1)*D;    
% pp.PoCHallMan  = Pc3D_Hall(r1,v1Man,C1,r2,v2,C2,HBR);


%% Output
if strcmpi(pp.orbit,'leo'); pp.figs.lla = false; end
pp.sqrMahaLim    = pp.majorIter(end).sqrMahaLim;
pp.sqrMahaMan    = smd;
pp.ipcMan        = ipc_tot;
pp.ipcManSingle  = ipc;
pp.valErr        = max(normOfVec(pp.validationAbsTraj(1:3,:) - ...
                                 pp.majorIter(end).absP))*pp.scaling(1);
%% Check validity of gradient SMD constraint
if pp.enableSmdGradConstraint
    pp = smdGradValidation(P,pp.ipcMan,smdLim,pp);                                 %
    for i = 1:N
        pp.smdGrad(:,i) = smd_grad(P(:,:,i),pp.validationTraj(1:3,i));         % [-] (3xN) Node-wise SMD gradient
    end
end
end