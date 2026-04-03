function [limsUp, limsLo] = limits(x,r,N,m,uMax,smdGradLim,gradOn,xi_dyn,xi_ca,nu_max,pp)
% Limits computes the vector of the limits for the optimization
% variables for the convex optimization in MOSEK.
%
% INPUT: xVec        = [-] (9,N) Expansion point for the state and control
%        N           = [-] (1,1) Number of nodes in the problem. 
%        m           = [-] (1,1) Number of variables per node
%        uMax        = [-] (1,1) Maximum control acceleration
%        smdGradLim  = [-] (N,1) Limit of the gradient of the SMD.
%        gradOn      = [bool] (N,1)  vector that activates the SMD gradient
%                      limits
%        xi_ca       = [-] (6,N) Non-linear parameter nu
%        nu_max      = [-] (1,1) Maximum allowed value for the NLI
%        pp          = [struct]  Postprocess structure
% 
% OUTPUT: limsUp = (mN,1) Upper limits of the optimization vector
%         limsLo = (mN,1) Lower limits of the optimization vector 
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
M                              = length(pp.secondary);
canFire                        = pp.canFire;                               % [bool] (N,1) Nodes in which the thrusters can fire
ind                            = pp.index;                                 % [struct]  Structure of the indeces
limsUp                         = 1*ones(m,N);                              % [-] (m,N) Initialize upper limits
limsLo                         = -1*ones(m,N);                             % [-] (m,N) Initialize lower limits
Dx                             = nu_max./xi_dyn(1:6,2:end);                % [-] (6,N) Maximum allowed state deviation from reference
Dr                             = nu_max./max(xi_ca,[],3);                  % [-] (3,M) or (3,N,M) Maximum allowed state deviation from reference
% Dr                             = nu_max./xi_ca;                     % [-] (3,M) or (3,N,M) Maximum allowed state deviation from reference
% Dr                             = 1./xi_ca;                     % [-] (3,M) or (3,N,M) Maximum allowed state deviation from reference
limsUp(ind.ctrlCone ,canFire)  = uMax;                                     % [-] (1,N) Maximum trhust
limsUp(ind.ctrlCone ,~canFire) = 0;                                        
limsUp(ind.ctrl     ,canFire)  = uMax;                   
limsUp(ind.ctrl     ,~canFire) = 0;
limsLo(ind.ctrlCone ,:)        = 0;                                        % Minimum thrust
limsLo(ind.ctrl     ,canFire)  = -uMax*~pp.justInTime;
limsLo(ind.ctrl     ,~canFire) = 0;
limsUp(ind.state    ,1)        = x(:,1);                                   % Fixed initial state
limsLo(ind.state    ,1)        = x(:,1); 
limsUp(ind.state  ,2:end)      = x(:,2:end) + Dx;                               % (6,N) Upper limits for the state
limsLo(ind.state  ,2:end)      = x(:,2:end) - Dx;                               % (6,N) Lower limits for the state  
if pp.pocConstr
    if pp.fastEncounter
        for j = 1:M
            i = pp.NCA(j);
            limsUp(ind.state(1:3),i) = min(limsUp(ind.state(1:3),i), ...
                                            x(1:3,i) + Dr(:,j));
            limsLo(ind.state(1:3),i) = max(limsLo(ind.state(1:3),i), ...
                                            x(1:3,i) - Dr(:,j));
        end
    end
end

limsUp(ind.vcCone ,:)     = 1e4;                                           % (1,N) Upper limits for the virtual control cone
limsLo(ind.vcCone ,:)     = 0;                                             % (1,N) Lower limits for the virtual control cone
limsUp(ind.vc ,:)         = 1e4;                                           % (*,N) Upper limits for the virtual controls and virtual buffers                                       % (1xN) Upper limits for the virtual control cone
limsLo(ind.vc ,:)         = -1e4;                                          % (*,N) Lower limits for the virtual controls and virtual buffers
if pp.enableSmdGradConstraint
    limsUp(ind.gradCone ,:) = 1e4;                                         % (1,N) Upper limit for the SMD gradient cone
    limsLo(ind.gradCone ,:) = 0;                                           % (1,N) Lower limit for the SMD gradient cone
    limsUp(ind.grad     ,:) = 1e4';                                        % (3,N) Upper limit for the SMD gradient
    limsLo(ind.grad     ,:) = -1e4';                                       % (3,N) Upper limit for the SMD gradient
end
if pp.smdSoft
    limsLo(ind.gradSlack ,:) = 0;                                          % (6,N) Lower limits for the SMD gradient slack variables
end
if pp.stationKeeping && pp.skSoft
    limsUp(ind.skSlack ,:) = 2*pi;                                         % (4,N) Upper limits for the SK slack variables
    limsLo(ind.skSlack ,:) = 0;                                            % (4,N) Lower limits for the SK slack variables
end 
if pp.enableHomotopy
    limsUp(ind.homotopy ,:) = uMax;                                        % (1,N) Upper limits for the homotopy variables
    limsLo(ind.homotopy ,:) = 0;                                           % (1,N) Lower limits for the homotopy variables 
end
limsUp = reshape(limsUp,m*N,1);                                            % (mN,1) reshape limits in a vector
limsLo = reshape(limsLo,m*N,1); 
if pp.enableSkTarget
    limsUp(ind.targetSlack) = 1e6;                                         % Add final target slack limits at the tail of the limit vectors
    limsLo(ind.targetSlack) = 0;
end  
end