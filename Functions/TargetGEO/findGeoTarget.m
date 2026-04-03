function [x0New,pp] = findGeoTarget(x0,tf,dt,nomLon,pp)
% FinGeoTarget computes the optimal target state for the end of the CA
% maneuver to stay inside the SK box for the required period of time. This
% is done using a convex problem that minimizes the violations of the SK
% box and the deviation from the original state.
% 
% INPUT:  x0     = [-] (6,1) Original target state
%         tf     = [-] (1,1) Final time of the propagation in which the
%                            violations need to be minimized
%         dt     = [-] (1,1) Time step
%         nomLon = [-] (1,1) Nominal longitude around which the SK box is
%                            defined
%         pp     = [struct] Postprocess structure
% 
% OUTPUT: x0New  = [-] (6,1) Optimized target state
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

%% Initialize variables
et          = pp.et + pp.N_forw*pp.dt;                                      % [-] (1,1) Ephemeris time of targeting point
N           = length(0:dt:tf);                                             % [-] (1,1) Number of nodes
err         = 1;                                                           % [-] (1,1) Initialize error
iter        = 0;                                                           % [-] (1,1) Initialize iteration counter
ind.state   = 1:6;                                                         % [-] (1,6) state indeces
ind.pos     = ind.state(1:3);                                              % [-] (1,3) position indeces
ind.skSlack = ind.state(end)   + (1:4);                                    % [-] (1,4) geodetic slack variables indeces
ind.vc      = ind.skSlack(end) + (1:12);                                   % [-] (1,12) virtual control indeces
m           = ind.vc(end);
ind.dx0     = m*N + (1:12);                                                   % [-] (1,3) position indeces
xi_max      = 1e-3;                                                        % [-] (1,1) Maximum allowed NLI, we want to stay as close as possible to the original point, so it is low
tol         = 1e-6;                                                        % [-] (1,1) Error tolerance
maxIter     = 10;                                                          % [-] (1,1) Maximum allowed number of iterations
newTraj     = nan;                                                         % Placeholder for the firsst iteration
%% Run iterations
while err > tol && iter < maxIter 
    iter = iter + 1;
    [dynMaps,llMaps,state,ll_const,nu] = ...
                        propagateTarget(N,x0,newTraj,dt,et,iter,pp);       % Perform the propagation
    pp.timeSubtr = pp.timeSubtr;
    if iter == 1
        xExp = state;                                                      % Expansion point for the first iteration taken from the forward propagation
    else
        xExp = xNew(1:6,:);                                                % Expansion point for iteration after the first one
    end
    xNew    = geoAltOpt(N,dynMaps,llMaps,nomLon,state,xExp, ...
                             ll_const,x0,xi_max,nu,ind,pp);                 % Perform the optimization
    err     = norm(xNew(1:6,1) - state(:,1));                              % [-] (1,1) Error computed on the initial state
    er((iter)) = err;
    newTraj = reshape(xNew(1:6,:),6*N,1);                                  % Reshape the new optimized trajectory
    slErr   = max(normOfVec(xNew(ind.skSlack,:)));                         % [-] (1,1) Maximum value of the geodetic slack variables
    
    % increase the weight of the minimization of the trust region if the 
    % slack variables are low enough. This allows to obtained a highly 
    % reliable solution once a proper point has been found
    if slErr < 1e-5                                                        
        xi_max = xi_max/10;
    end
end
x0New = xNew(1:6,1);                                                       % [-] (6,1) Optimized target state
end