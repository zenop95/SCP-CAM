function xNew = geoAltOpt(N,dynMaps,llMaps,nomLon,state,xExp,ll_const,x0, ...
                                                xi_max,nu,ind,pp)
% geoAltOpt solves the convex problem to find the initial position and 
% velocity that minimize the deviations in latitude and longitude for a 
% GEO satellite.
%
% INPUT:  dynMaps   = [-] (6,6,N) Linear maps for the dynamics
%         llMaps    = [-] (2,3,N) Linear maps from ECI to geodetic   
%         nomLon    = [-] (1,1)   Nominal value of the longitude
%         state     = [-] (6,N)   constant part of the propagated state
%         x_exp     = [-] (6,N)   Expansion point for the linearization
%         ll_const  = [-] (2,N)   Constant part of the geodetic coordinates
%         N         = [-] (1,1)   Number of nodes in the optimization
%         xi_max    = [-] (1,1)   Maximum allowed value for the NLI
%         w         = [-] (1,1)   Weight of the NLI cone in cost function
%         nu        = [-] (6,N)   Nonlinearity parameter nu
%         ind       = [struct]    Structure of the indeces for each node
%         pp        = [struct]    Postprocess structure
% 
% OUTPUT: xNew      = [-] (6,N) Optimized vector
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
m     = ind.vc(end);                                                       % [-] (1x1) Number of variables per node
nOpt  = m*N+12;                                                               % [-] (1x1) Total number of optimization variables
x_exp = [reshape([xExp; zeros(m-6,N)],m*N,1);zeros(12,1)];                              % [-] (mOptx1) Reshaped expansion point
%% Objective function
c                  =  zeros(m,N);                                     % Initialize the coefficients of the objective function
c(ind.skSlack , :) = 10;                                                   % Coefficient of the geodetic slack variable to minimize violations
c(ind.vc      , :) = 1e3;                                                  % Coefficient of the virtual control
c                  = [reshape(c,m*N,1); 1*ones(12,1)];                                    % Reshape into vector

%% Limits
limsUp                 = nan(m,N);                                         % Initialize the upper limits of the optimization variables
limsLo                 = nan(m,N);                                         % Initialize the lower limits of the optimization variables
Dx                     = xi_max./nu;                                       % Maximum allowed deviations to the state variables
limsUp(ind.state   ,:) = xExp + Dx;                                        % Upper limits of the state variables
limsLo(ind.state   ,:) = xExp - Dx;                                        % Lower limits of the state variables
limsUp(ind.skSlack ,:) = 2*pi;                                             % Upper limits of the slack variables
limsLo(ind.skSlack ,:) = 0;                                                % Lower limits of the slack variables
limsUp(ind.vc      ,:) = 1;                                                % Upper limits of the virtual control variables
limsLo(ind.vc      ,:) = 0;                                                % Lower limits of the virtual control variables

%% equality: continuity of the trajectory
A_dyn  = zeros(6*(N-1), nOpt);
A_dyn1 = A_dyn;
for i = 2:N
    A_dyn(1+6*(i-2):6*(i-1),  ind.state  + m*(i-2)) = dynMaps(:,:,i);      % Coefficients of the state variables
    A_dyn(1+6*(i-2):6*(i-1),  ind.vc     + m*(i-2)) = [eye(6), -eye(6)];   % Coefficients of the virtual control variables
    A_dyn(1+6*(i-2):6*(i-1),  ind.state  + m*(i-1)) = -eye(6);             % Coefficients of the propagated state variables
    A_dyn1(1+6*(i-2):6*(i-1), ind.state  + m*(i-2)) = dynMaps(:,:,i);      % Matrix to perform the propagation of the expansion point
end
A_dyn      = sparse(A_dyn);                                                % Sparse matrix of the coefficients of the dynamics constraint
A_dyn1     = sparse(A_dyn1);
const_part = reshape(state(:,2:end),6*(N-1),1);                            % Constant part of the propagation
b_dyn      = A_dyn1*x_exp - const_part;                                    % Residual of the linearization

%% SK constraint
A_sk      = zeros(2*N,nOpt);
expLatLon = nan(2,N);           
for i = 1:N
    A_sk(1+(i-1)*2:2*i , ind.pos     + (i-1)*m) = llMaps(:,:,i);           % Coefficients of the ECI position variables
    A_sk(1+(i-1)*2:2*i , ind.skSlack + (i-1)*m) = [1 -1 0 0; 0 0 1 -1];    % Coefficients of the slack variables
    expLatLon(:,i) = llMaps(1:2,:,i)*xExp(1:3,i);                          % Conversion of the expansion point into geodetic coordinates
end
A_sk = sparse(A_sk);                                                       % Sparse matrix of the coefficients of the SK constraint
b_sk_up = [deg2rad(pp.skDev(1)); 
        nomLon + deg2rad(pp.skDev(2))].*ones(2,N) + expLatLon - ll_const;  % Residual of the linearization + upper limit
b_sk_lo = [-deg2rad(pp.skDev(1)); 
        nomLon - deg2rad(pp.skDev(2))].*ones(2,N) + expLatLon - ll_const;  % Residual of the linearization + lower limit
%% Dx0 constraint
A_dx  = zeros(6,nOpt);
A_dx(:, ind.dx0(1:6))    = eye(6);           
A_dx(:, ind.dx0(7:12))   = -eye(6);           
A_dx(:, ind.state)  = -eye(6);   
A_dx = sparse(A_dx);                                                       % Sparse matrix of the coefficients of the constraint
b_dx = x0;  % upper-lower limit

%% 
prob              = struct();
prob.c            = c;                                                               
prob.a            = [A_dyn; A_sk; A_dx];
prob.buc          = [b_dyn; reshape(b_sk_up,2*N,1); b_dx];
prob.blc          = [b_dyn; reshape(b_sk_lo,2*N,1); b_dx];
prob.bux          = [reshape(limsUp,m*N,1); 10*ones(12,1)];
prob.blx          = [reshape(limsLo,m*N,1); zeros(12,1)];

%% Solve the problem
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-11;
param.MSK_DPAR_INTPNT_TOL_PFEAS = 1e-11;
[~,res] = mosekopt('minimize echo(0)',prob,param);
% param.MSK_IPAR_INFEAS_REPORT_AUTO = 1;
% [~,res] = mosekopt('minimize',prob,param);
xNew  = res.sol.itr.xx;
xNew  = reshape(xNew(1:end-12),m,N);
end