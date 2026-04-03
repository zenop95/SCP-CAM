function nor = normOfVec(vec)
%normOfVec accepts a succession of column vectors and computes the norm 
% over each column. The ouput is a row vector of the norms.

% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

N = size(vec,2);
nor = nan(1,N);
for i = 1:N
    nor(i) = norm(vec(:,i));
end
end

