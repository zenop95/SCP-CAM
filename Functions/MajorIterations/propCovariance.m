function [P,pp,C_p] = propCovariance(stm,cartMaps,N,sum,pp)
% PropCovariance propagates covariancefor the duration of the propagation
% window. The propagation must be performed with STMs in cartesian
% coordinates
% 
% INPUT: STM         = [-]    (6,6,N) state trasition matrices
%        cartMaps    = [-]    (6,6,N) maps that transform state to cartesian
%        N           = [-]    (1,1)   Number on nodes in the propagation
%        sum         = [bool] (1,1)   If true, the secondary covariance is 
%                                     summed to the pimray, if false it is 
%                                     not. It is kept false for the back 
%                                     propagation of the covariance.
%        pp          = [struct]       Postprocess structure.
%        
% OUTPUT: C          = [-] (6,6,N) State covariance for each node.
%         P          = [-] (3,3,N) Position covariance for each node.
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
M   = length(pp.secondary);
C   = nan(6,6,N,M);
C_p = nan(6,6,N);
cartStm = nan(6,6,N);
if ~strcmp(pp.obj,'miss_distance')
    cartStm(:,:,1) = eye(6);
    for i = 2:N
        cartStm(:,:,i) = cartMaps(:,:,i)*stm(:,:,i)/cartMaps(:,:,i-1); % [-] (6,6) STM in cartesian coordinates
    end
    if ~isfield(pp,"majorIter")
        C_p(:,:,pp.N_back+1) = pp.primary.C0;
        for i = pp.N_back+1:-1:2
            stm_inv    = inv(cartStm(:,:,i));
            C_p(:,:,i-1) = stm_inv*C_p(:,:,i)*stm_inv';
        end
        for i = pp.N_back+2:N
            C_p(:,:,i) = cartStm(:,:,i)*C_p(:,:,i-1)*cartStm(:,:,i)';      % [-] (6,6) Propagate node i-1 to obtain covariance in node i
        end
        pp.primary.C0 = C_p(:,:,1);
    else
        C_p(:,:,1) = pp.primary.C0;
        for i = 2:N
            C_p(:,:,i) = cartStm(:,:,i)*C_p(:,:,i-1)*cartStm(:,:,i)';      % [-] (6,6) Propagate node i-1 to obtain covariance in node i
        end
    end
end
for j = 1:M
    if strcmp(pp.obj,'miss_distance')
        C(:,:,:,j) = repmat(eye(6),1,1,N);                                 % [-] (6,6,N) If the miss distance is used, the covariance is set to eye to obtain the squared of the distance in the SMD formula
    else
        if sum
            C(:,:,:,j) = C_p + pp.secondary(j).covariance;                 % [-] (6,6,N) Add covariance of the secondary s/c (Gaussian)
        else 
            C(:,:,:,j) = C_p;
        end
    end
end
P = C(1:3,1:3,:,:);                                                          % [-] (3,3,N) Extract position covariance
end