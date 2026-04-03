function [gradOn,gradLim,pp] = smdGradParams(iter,gradLim,smdLim,smd,ipc,pp)
% SmdGradParams constructs the parameters for the SMD gradient constraint
%
% INPUT: iter     = [-] (1,1) Current iteration counter
%        gradLim  = [-] (N,1) Initialized array for SMD gradient limit
%        ipc      = [-] (N,1) Node-wise value of the IPC
%        pp       = [struct] Postprocess structure.
% 
% OUTPUT: gradOn  = [bool] (N,1) Toggles the SMD gradient constraint
%         gradLim = [-]    (N,1) Array of SMD gradient limit
%         pp      = [struct] Postprocess structure.
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
% The difference in Pic at a distance equal to HBR must be lower than
% tol*IPCmax. This difference is computed multiplying the norm of the
% gradient of the SMD by the Hard Body Radius
if pp.gradLimContinuation
    tol = pp.maxIpcDeviation*(10-9*exp(iter-7)./(exp(iter-7)+0.05)); % Continuation function (sigmoid)
%     tol = pp.maxIpcDeviation+(1-iter/(pp.iterMaxMaj)); % Continuation function (sigmoid)
else
    tol = pp.maxIpcDeviation;                                              % [-] (1,1) Required tolerance on the IPC variation at a distance equal to HBR from nominal relative position
end
gradOn = findLocalMax(ipc);
% gradOn = zeros(pp.N,1);                                                  % [bool] (N,1) Initialize variable that toggles the SMD gradient constraint
gradOn(ipc < pp.ipcLim(1)*(1-pp.sdmGradTrig)) = 0;                         % [bool] (N,1) Toggles the SMD gradient constraint
for i = pp.NCA0:pp.NCAf
%         gradOn(i) = ipc(i) == max(ipc(pp.NCA0:pp.NCAf)); % activate only on maximum IPC point
%     gradOn(i) = ipc(i) >= pp.ipcLim(1)*(1-pp.sdmGradTrig);                 % [bool] (N,1) Toggles the SMD gradient constraint
%     gradOn(i) = smd(i) <= smdLim(i)*(1+pp.sdmGradTrig);
    if gradOn(i)
        gradLim(i) = 2*tol/pp.secondary.HBR*pp.ipcLim/ipc(i);                        % [-] (N,1) Array of SMD gradient limit
    end
end
if any(gradOn)
    pp.gradScale = max(gradLim(gradOn==1));                                    % [-] (1,1) Maximum value of the SMD gradient to scale the constraint in the optimization
else
    pp.gradScale = 1;
end
end