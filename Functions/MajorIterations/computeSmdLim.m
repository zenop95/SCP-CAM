function smdLim = computeSmdLim(P,r,highRisk,pp)
% propCovariance propagates covariance and yields limits for IPC and SMD
%
% INPUT: P       = [-] (3,3,N) Covariance matrix for each node.
%        pp      = [struct]    Postprocess structure.
%         
% OUTPUT: ipcLim = [-] (N,1) Limit IPC for each node
%         smdLim = [-] (N,1) Limit SMD for each node
%
% Documentation: https://docs.mosek.com/9.2/toolbox/tutorial-cqo-shared.html 
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
N   = pp.N;
obj = pp.obj;
M   = length(pp.secondary);
if strcmp(obj,'miss_distance')
    if pp.enableBplaneAvoidance
        smdLim = pp.dMin^2*ones(1,M);                                      % [-] (N,1) In the case of the miss distance constraint, the limit is the square of the required miss distance
    else
        smdLim = pp.dMin^2*ones(N,M);                                      % [-] (N,1) In the case of the miss distance constraint, the limit is the square of the required miss distance
    end
else
    for j = 1:M
        lim = pp.lims(j)/pp.w(j);
        HBR    = pp.secondary(j).HBR;
        if pp.enableBplaneAvoidance
            i         = pp.NCA(j);
            e2b2      = pp.e2b([1 3],:,j);     
            Pb        = e2b2*P(:,:,i,j)*e2b2';                             % [-] (2,2) 2D covariance in the B-plane
            detP      = det(Pb);                                           % [-] (1,1) determinant of the 2D covariance in the B-plane
            if strcmpi(pp.PoCType,'maximum')
                smdLim(j) = HBR^2/(exp(1)*sqrt(detP)*lim);                 % [-] (1,1) SMD limit computed with Maximum Alfriend and Akella's formula applied to PC
            elseif strcmpi(pp.PoCType,'chan')
                smdLim(j)= PoC2SMD(Pb, HBR, lim, 5, 1, 1e-3, 200);         % [-] (1,1) SMD limit computed with Chan's formula applied to PC
            else
                smdLim(j) = -2*log(2*lim*sqrt(detP)/HBR^2);                % [-] (1,1) SMD limit computed with Alfriend and Akella's formula applied to PC
            end
        else
            detP = nan(N,1);
            for i = 1:N
                detP(i)  = det(P(:,:,i,j));                                % [-] (N,1) Node-wise determinant of the positional covariance matrix
            end
            switch lower(pp.ipc_type)
                case 'constant'
                    smdLim(:,j) = -2*log(3*lim*sqrt(pi*detP)/ ...
                                                        (sqrt(2)*HBR^3));  % [-] (N,1) SMD limit computed with Alfriend and Akella's formula applied to IPC
                case 'cuboid'
                    for i = 1:N
                        if highRisk(i)
                            smdLim(i,j) = ipcCuboid2smd(P(:,:,i,j),r(:,i),HBR,lim,1e-5);     % [-] (N,1) SMD limit computed with Alfriend and Akella's formula applied to IPC
                        end
                    end
                case 'maximum'
                    smdLim(:,j) = (sqrt(2)*HBR)^3./(3*exp(1)* ...
                                                sqrt(pi*detP)*lim);        % [-] (N,1) SMD limit computed with the maximum risk formula
                otherwise
                    error('Invalid IPC type.')
            end
        end
        % If the IPClim is too high, the d_m^2 computed with these formulas can
        % be negative. Since these nodes do not influence the maneuver, the SMD
        % limit is set as very low in them.
        smdLim(smdLim <= 0) = 1e-5;                                        % [-] (N,1) 
    end
end
end