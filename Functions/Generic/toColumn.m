function col = toColumn(vec)
%TOCOLUMN accepts a vector in input and outputs the same vector as a column
%vector.

% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

validateattributes(vec,"double","vector")
n   = length(vec);
col = reshape(vec,[n 1]);
end

