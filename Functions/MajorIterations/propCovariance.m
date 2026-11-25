function [P,pp,C_p] = propCovariance(stm,cartMaps,N,combineSecondary,pp)
% propCovariance propagates covariance for the duration of the propagation
% window. The propagation must be performed with STMs in cartesian
% coordinates
% 
% INPUT: STM         = [-]    (6,6,N) state trasition matrices
%        cartMaps    = [-]    (6,6,N) maps that transform state to cartesian
%        N           = [-]    (1,1)   Number on nodes in the propagation
%        combineSecondary  = [bool] (1,1)   If true, the secondary covariance is 
%                                     summed to the primary, if false it is 
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
Np = numel(pp.primary);
Ns = numel(pp.secondary);
K  = Np*Ns;

C_p = nan(6,6,N,Np);
P   = zeros(3,3,N,K);

% Propagate primary covariance for every mixand
if ~strcmp(pp.obj,'miss_distance')
    for a = 1:Np
        cartStm = zeros(6,6,N);
        cartStm(:,:,1) = eye(6);

        % Convert each primary mixand STM to Cartesian coordinates
        for i = 2:N
            cartStm(:,:,i) = cartMaps(:,:,i,a) * stm(:,:,i,a) / cartMaps(:,:,i-1,a);
        end

        if ~isfield(pp,"majorIter")
            % Initial covariance is defined at TCA
            C_p(:,:,pp.N_back+1,a) = pp.primary(a).C0;

            % Backward propagation from TCA to node 1
            for i = pp.N_back+1:-1:2
                stm_inv = inv(cartStm(:,:,i));

                C_p(:,:,i-1,a) = stm_inv*C_p(:,:,i,a)*stm_inv';
            end

            % Forward propagation from TCA to node N
            for i = pp.N_back+2:N
                C_p(:,:,i,a) = cartStm(:,:,i)*C_p(:,:,i-1,a)*cartStm(:,:,i)';
            end

            % Store covariance at the start of the window
            pp.primary(a).C0 = C_p(:,:,1,a);

        else
            % Subsequent iterations: covariance is known at node 1
            C_p(:,:,1,a) = pp.primary(a).C0;
            for i = 2:N
                C_p(:,:,i,a) = cartStm(:,:,i)*C_p(:,:,i-1,a)*cartStm(:,:,i)';
            end
        end
    end
end


for a=1:Np
    for c=1:Ns
        k=(a-1)*Ns+c;
        if combineSecondary
            Crel=C_p(:,:,:,a)+pp.secondary(c).covariance;
        else
            Crel=C_p(:,:,:,a);
        end
        P(:,:,:,k)=Crel(1:3,1:3,:);
        % P(:,:,:,k) = 0.5*(Pk + permute(Pk,[2 1 3]));
    end
end
 
end