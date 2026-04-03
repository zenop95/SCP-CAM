function [sqrMaha, highRisk, ipc] = checkMahaCond(relTraj,P,sqrMahaLim,pp)
N = pp.N;
HBR = repmat(pp.HBR,N,1);
ipc_type = repmat(string(pp.ipc_type),N,1);
[sqrMaha, highRisk, ipc] = arrayfun(@(relTraj,P,sqrMahaLim,HBR,ipc_type) ...
                    function_handle(relTraj,P,sqrMahaLim,HBR,ipc_type),...
                    relTraj,P,sqrMahaLim,HBR,ipc_type);
end

function [sqrMaha, highRisk, ipc] = function_handle(relTraj,P,sqrMahaLim,HBR,ipc_type)
    r = relTraj(1:3);
    P = reshape(P,3,3);
    sqrMaha = dot(r,P\r);
    highRisk(sqrMaha < sqrMahaLim) = 1;
    switch obj
        case 'max_risk'
            ipc = maximumIpc(r,P,HBR);
        case 'risk'
            switch lower(ipc_type)
                case 'constant'
                ipc = constantIpc(r,P,HBR);
                case 'cuboid'
                ipc = cuboidIpc(r,P,HBR);
                case 'numeric'
                ipc = numericIpc(r,P,HBR);
            end
    end
end
