function pp = smdGradValidation(P,ipc,smdLim,pp)
%UNTITLED Summary of this function goes here
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

N = pp.N;
smdGrad           = nan(3,N);
pp.deltaIpc       = nan(N,1);
smdGradNormalized = nan(3,N);
for i = 1:N
    [V,D] = eig(P(:,:,i));
    Pc    = V'*P(:,:,i)*V;
    r     = V*pp.validationTraj(1:3,i);
    smdGrad(:,i)           = 2*r./diag(Pc);
    smdGradNormalized(:,i) = V'*normalize(smdGrad(:,i),'norm');
end
deviation = pp.validationTraj(1:3,:) - smdGradNormalized*pp.secondary(1).HBR; 
pp.ipcAtOffset = computeIpc(deviation,P,pp);
% Now compute the same thing using the linear formula DP = P*|grad|*HBR/2
for i = 1:N
    if pp.enableBplaneAvoidance
        pp.deltaIpc(i) = ipc*norm(smdGrad)*pp.secondary(1).HBR/2;
    else
        pp.deltaIpcLin(i) = ipc(i)*norm(smdGrad(:,i))*pp.secondary(1).HBR/2;
    end
end
pp.deltaIpc = pp.ipcAtOffset - ipc;
pp.ipcLin = pp.deltaIpc + ipc;
end