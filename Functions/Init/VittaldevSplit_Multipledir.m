function [xi, Pi, wi] = VittaldevSplit_Multipledir(x0, P0, dir, N)
%
% Multidirectional Vittaldev Gaussian Mixture Splitting
%
% DESCRIPTION:
% Recursively splits an initial Gaussian distribution along multiple
% specified directions using the univariate splitting library.
%
% PROTOTYPE:
%   [xi, Pi, wi] = VittaldevSplit(x0, P0, dir, N)
%
% INPUTS:
%   x0    [6x1]   Initial state vector
%   P0    [6x6]   Initial covariance matrix
%   dir   [6xm]   Splitting directions
%   N     [1xm]   Number of components at each splitting stage
%
% OUTPUTS:
%   xi    [6xK]   Gaussian component means
%   Pi    [6x6xK] Gaussian component covariance matrices
%   wi    [Kx1]   Gaussian component weights
%
% where K = prod(N).
%
%% Input processing

x0 = x0(:);
N = N(:).';

if size(dir,1) ~= length(x0)
    error('Splitting directions must have six rows.');
end

if size(dir,2) ~= numel(N)
    error('One splitting count is required per direction.');
end

if any(N < 1) || any(N ~= round(N))
    error('Splitting counts must be positive integers.');
end

%% Initial Gaussian

xi = x0;
Pi = (P0 + P0')/2;
wi = 1;

%% Recursive multidirectional splitting

for k = 1:numel(N)

    Nk = N(k);
    dk = dir(:,k);

    if norm(dk) == 0
        error('Splitting direction cannot be zero.');
    end

    dk = dk/norm(dk);

    % Retrieve univariate splitting locations and weights
    [mu_i, w_i] = SplitLib(Nk);

    mu_i = mu_i(:);
    w_i = w_i(:);

    if numel(mu_i) ~= Nk || numel(w_i) ~= Nk
        error('SplitLib returned an incorrect number of components.');
    end

    nParent = size(xi,2);
    nChild = nParent*Nk;

    % Allocate next generation
    xi_new = zeros(6,nChild);
    Pi_new = zeros(6,6,nChild);
    wi_new = zeros(nChild,1);

    % Split every Gaussian from the previous generation
    for j = 1:nParent

        x_parent = xi(:,j);
        P_parent = Pi(:,:,j);

        P_parent = (P_parent + P_parent')/2;

        % Symmetric square-root decomposition
        [V,D] = eig(P_parent);

        eigValues = diag(D);

        if any(eigValues <= 0)
            error('Parent covariance is not positive definite.');
        end

        S = V*diag(sqrt(eigValues))*V';

        % Direction-dependent displacement
        shift_mean = dk/norm(S\dk);

        % Child covariance
        P_child = P_parent + ...
            (1/Nk - 1)*(shift_mean*shift_mean');

        P_child = (P_child + P_child')/2;

        % Generate children
        for i = 1:Nk

            index = (j-1)*Nk + i;

            xi_new(:,index) = ...
                x_parent + mu_i(i)*shift_mean;

            Pi_new(:,:,index) = P_child;

            wi_new(index) = ...
                wi(j)*w_i(i);

        end
    end

    % Update Gaussian mixture
    xi = xi_new;
    Pi = Pi_new;
    wi = wi_new;

end

%% Weight normalization

wi = wi/sum(wi);

end
