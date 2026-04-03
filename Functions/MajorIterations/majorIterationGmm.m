function [newTraj,majorIter,xi_max,breakFlag,pp] = ...
                              majorIterationGmm(oldTraj,xi_max,iter,k,G_db,pp)
% MajorIteration runs the major iteration process for the update of the
% dynamics of the solution.
%
% INPUT: oldTraj = [-]    (6,N) Relative trajectory from previous iteration
%        iter    = [-]    (1,1) Itaration counter
%        hTrue   = [bool] (N,1) 1 if homotopy is respected, 0 otherwise
%        pp      = [struct]     Postprocess structure
%        
% OUTPUT: newTraj   = [-]    (6,N) Optimized trajectory
%         majorIter = [struct]     Structure containing output major iteration
%         breakFlag = [bool] (1,1) Flag that stops the simulation
%         pp        = [struct]     Postprocess structure
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

N       = pp.N;
newTraj = nan;
if iter == 1 && pp.N_back > 0 
    [pp.x_s, pp.C0p] = propBack(pp);                                       % [-] (6,1) (6,6) Back propagate from conjunction if required
end

%% GMM split
[V,D] = eig(pp.C0p);
[a,b] = sort(diag(D),'descend');
D = diag(a);
V = V(:,b);
C0c   = D;
xc    = V'*pp.x_s;
dir   = [1; zeros(5,1)];
n     = pp.n;
[xic, Pic, w] = VittaldevSplit(xc, C0c, dir, n);
for j = 1:n
    xi(:,j) = V*xic(:,j);
    Pi(:,:,j) = V*Pic(:,:,j)*V';
end

GmmMean = nan(6,N,n);
wMean   = nan(6,N,n);
gmmPos   = nan(6,N,n);
GmmCov  = nan(3,3,N,n^2);
wCov    = nan(3,3,N,n^2);
GmmRelTraj = nan(6,N,n^2);
l = 0;
A = nan(6,9,N,n);
G = nan(3,6,N,n);
ll = nan(3,N,n);
for j = 1:n
    %% Propagate and extract maps and constant parts 
    [dynMaps,STM,llMaps,cartMaps,state,cart,lla,nu,pp] = ...
                                        propagateDA(iter,oldTraj,false,pp);
%     pp.nomLon = -2.7066566;
    A(:,:,:,j) = dynMaps;
    G(:,:,:,j) = llMaps;
    ll(:,:,j) = lla;
    pp.nomLon = lla(2,1);                                                % [-] (1,1) Nominal value of longitude for GEO station keeping
    relTraj = cart - pp.secondary.cart;                                    % [-] (6,N) Relative trajectory
    uP = zeros(N,1);                                                       % Kept but not used, if needed use relStmProp
    GmmMean(:,:,j) = state;
    wMean(:,:,j)   = state*w(j);
    %% Covariance propagation
    for k = 1:n
        l = l+1;
        if pp.propCovariance
            [~,P] = propCovariance(STM,cartMaps,pp.N,true,k,pp);                     % [-] (6,6,N) Propagated positional covariance matrix
        else
            P = nan(3,3,N);
        end
        GmmCov(:,:,:,l)   = P;
        wCov(:,:,:,l)     = P*w(j)*w(k);
        GmmRelTraj(:,:,l) = GmmMean(:,:,j)-pp.secondary.gmmMean(:,:,k);
        weights(l) = w(j)*w(k);
    end
end
cart    = sum(wMean,3);
dynMaps = mean(A,4);
llMaps  = mean(G,4);
lla     = mean(ll,3);
for j = 1:n
    gmmPos(:,:,j) = cart*(1-w(j)) + GmmMean(:,:,j)*w(j);
end
%% Check where SMD condition is violated
smdLim             = computeSmdLimGmm(GmmCov,weights,pp);                          % [-] (N,1)   Compute the node-wise SMD limits
[ipc,smd,highRisk] = computeIpcGmm(GmmRelTraj,GmmCov,weights,smdLim,pp);           % [-] (N,1)x3 Compute node-wise SMD and IPC
ipc = sum(ipc,1);
%% GEO target state
if pp.enableSkTarget && iter == 1
    if strcmpi(pp.orbit,'geo') && ~pp.loadTarget
        dt     = 5000/pp.Tsc;                                              % [-] (1,1) Time step for the computation of the target state
        tf     = 14*pp.T;                                                  % [-] (1,1) Ending time of propagation for the computation of the target state
        cart_f = pp.x2cart(state(:,end));                                  % [-] (6,1) Final cartesian state unmaneuvered
        x0New  = pp.cart2x(findGeoTarget(cart_f,tf,dt,pp.nomLon,pp));      % [-] (6,1) Target state
    elseif strcmpi(pp.orbit,'geo') && pp.loadTarget
        x0New = pp.cart2x(load('skTarget.mat').x0New);                     % [-] (6,1) Load a previously defined target point, stored in the file 'skTarget.dat'
    else
        x0New = state(1:6,end);                                            % [-] (6,1) In the LEO case the final target is the return to the nominal orbit
    end
    pp.skTarget   = x0New;
    pp.cartTarget = pp.x2cart(x0New);                                      % [-] (6,1) Final cartesian target
end
%% SMD gradient limit constraint
gradLim = 1e5*ones(N,1);
gradOn = zeros(N,1);
pp.gradScale = 1;
if pp.enableSmdGradConstraint && iter > 1                                  % First iteration with no constraint
    [gradOn,gradLim,pp] = smdGradParams(iter,gradLim,smdLim,smd,ipc,pp);   % [-] (N,1)x2 Compute the node-wise limit and set the nodes in which to activate it
end
pp.gradOn = gradOn;

%% Ballistic things
if iter == 1
    pp.ballisticTraj    = cart;                                            % [-] (6,N) Cartesian trajectory of the primary
    pp.ballisticRelTraj = cart - pp.secondary.cart;                        % [-] (6,N) Cartesian relative trajectory
    pp.ballisticState   = state;                                           % [-] (6,N) Node-wise state of the primary
    pp.ballisticIpc     = ipc;                                             % [-] (N,1) Node-wise IPC
    pp.ballisticSmd     = smd;                                             % [-] (N,1) Node-wise SMD
    for i = 1:N
        pp.ballisticSmdGrad(:,i) = smd_grad(P(:,:,i),relTraj(1:3,i));      % [-] (3,N) Node-wise SMD gradient
    end
end
%% Solve the problem using successive convexification (minor iterations)
constPos = cart(1:3,:);                                                    % [-] (3,N) Constant part of the cartesian position
refute = true;
while refute
    if iter > 1
        uP = pp.majorIter(end).uP;
%         [xi_max,J,refute] = updateXiMax(pp,state,oldTraj,iter);
    end
    [xNew,majorIter,breakFlag] = linConvexSolve(dynMaps,llMaps,cartMaps, ...
       state,k,G_db,uP,lla,constPos,GmmCov,gmmPos,smdLim,highRisk,gradLim,gradOn,nu,xi_max,pp);
    newTraj     = xNew(1:6,:);
%     if iter == 1
        refute = false; 
        J = sum(majorIter.minorIter(end).uNorm)*pp.ctrlWeight;
%     end
    if breakFlag; warning('Infeasible solution'); return; end
end
majorIter = majBuild(majorIter,lla,P,smdLim,STM,nu,xi_max,J,gradLim,pp);           % [struct] Build the structure of the major iteration
end