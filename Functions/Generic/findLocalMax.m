function res = findLocalMax(x)
%findLocalMax yeilds a vector in which the non-zero entries correspond to
%the local maxima of the input
% INPUT: x = [] Nx1 vector
% 
% OUTPUT: res = [] Nx1 vector with one entries in positions in which x has
% a local maximum
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
N = length(x);
res = zeros(N,1);
res(1) = x(1) > x(2); res(end) = x(end) > x(end-1);
for i = 2:N-1
  res(i) = x(i) > x(i-1) && x(i) > x(i+1);
end
end