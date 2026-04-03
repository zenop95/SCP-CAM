function pp = defIndConv(pp)
% DefIndConv defines the indeces of the variables per each node to set
% the SOCP in MOSEK.
%
% INPUT: pp  = [struct] structure with initialization parameters
%
% OUTPUT: pp = [struct] structure with initialization parameters
%
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
M = length(pp.secondary);
i.state    = 1:6;                                                          % state of the s/c
if pp.justInTime 
    i.ctrl = 7;                                                            % throttle (just-in-time CA)
else
    i.ctrl = 7:9;                                                          % control (full CA)
end
i.ctrlCone = i.ctrl(end) + 1;
m = i.ctrlCone(end);
i.grad = []; i.gradCone = []; i.gradVc = []; i.gradSlack = []; i.vcCone = [];
i.skVc = []; i.skSlack = []; i.homotopy = [];  i.stateVc = []; i.caVc = [];
i.altVc = [];
i.stateVc = m + (1:6);                                             % virtual control acting on the state
i.caVc    = i.stateVc(end) + (1:M);                                        % virtual buffer acting on the CA constraint
i.vcCone  = i.caVc(end) + 1;                                               % cone on the virtual control variables
m         = i.vcCone;                                                      % number of variables per node

if pp.enableSmdGradConstraint
    i.grad     = m + (1:3);                                                % SMD gradient
    i.gradCone = i.grad(end) + 1;                                          % cone on the SMD gradient
    i.gradVc = i.gradCone + 1;                                         % virtual buffer acting on the SMD gradient constraint
    m        = i.gradVc;
    if pp.smdSoft
        i.gradSlack = m + (1:6);                                           % slack acting on the SMD gradient for soft constraint
        m           = i.gradSlack(end);
    end
end
if pp.stationKeeping
    i.skVc = m + (1:2);                                                % virtual buffer acting on the station keeping constraint
    m      = i.skVc(end);
    if pp.skSoft
        i.skSlack = m + (1:4);                                             % slack acting on the station keeping constraint
        m         = i.skSlack(end);
    end
end
if pp.altSk
    i.altVc = m + 1;
    m       = i.altVc(end);
end
if pp.enableHomotopy
    i.homotopy = m + 1;                                                    % homotopy variable for minimum control action
    m          = i.homotopy;
end
if pp.enableSkTarget
    i.targetSlack = m*pp.N + (1:12);                                       % slack for soft constraint on final target 
end
i.vc = [i.stateVc,i.caVc,i.gradVc,i.skVc,i.altVc];
pp.index = i;                                                              % create final structure
pp.m     = m;
end