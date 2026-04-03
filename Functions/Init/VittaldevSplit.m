function [xi, Pi, wi] = VittaldevSplit(x0, P0, dir, N)
%
% Vittaldev Splitting
%
%DESCRIPTION:
%This code provides the splitting algorithm in order to create the Gaussian
%Mixtures.
%
%PROTOTYPE
%   [xi, Pi, wi] = VittaldevSplit(x0, P0, dir, N)
%
%--------------------------------------------------------------------------
% INPUTS:
%   x0         [6x1]       State Vector              [km], [km][s-1] (or [-])
%   P0         [6x6]       Covariance Matrix         [km2], [km2][s-2] (or [-])
%   dir        [6x1]       Direction of Splitting    [-]
%   N          [1x1]       Number of Mixtures        [-]
%--------------------------------------------------------------------------
% OUTPUTS:
%   xi          [Nx6]      Mixtures States            [km], [km][s-1] (or [-])
%   Pi          [6x6xN]    Mixtures Cov. Matrices     [km2], [km2][s-2] (or [-])
%   wi          [1xN]      Weights                    [-]
%--------------------------------------------------------------------------
%
%NOTES:
% - The Direction of Splitting can be defined with various tachniques:
%    * Spectral Decomposition (eigenvector associated to maximum value
%      aigenvalue)
%    * Square Root Matrix
% - The input "N" shall assume one of the following values (to meet the
%   "uni_mat.mat" map keying):
%   11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 3, 31, 33, 35, 37, 39, 5, 7, 9
% - The states and cov.mat. can also be gieven in Synodic.
%
%CALLED FUNCTIONS:
% SplitLib
%
%UPDATES:
% 2022/08/10, Luigi De Maria: added header and comments
%
%REFERENCES:
% (none)
%
%AUTHOR(s):
%unknown
%

%% Main Code

%Extraction of State and Splitting Direction
x0  = x0(:);
dir = dir(:);

%Eigenproblem of Covariance Matrix
[v,w] = eig(P0);   	%Eig decomposition of covariance
n     = length(x0);	%Problem dimension

sqrtD = sqrt(w);
S = v*(sqrtD*v');

astar_unnorm = S\dir;
shift_mean = dir / norm(astar_unnorm);

p_i = 1/N;
P_i = P0 + (p_i - 1.0) * (shift_mean*shift_mean');

%Load Maps for Splitting
[mu_i, w_i] = SplitLib(N);

%Memory Allocation
xi = zeros(n,N);   %States
Pi = zeros(n,n,N); %Covariance Matricies
wi = zeros(N,1);   %Weights

%Mixtures Data Creation
for i = 1 : N
    xi(:,i) = reshape(x0,[n,1]) + mu_i(i) * shift_mean;
    Pi(:,:,i) = P_i;
    wi(i) = w_i(i);
end

end