function plot3d(a,varargin)
%Plots in 3d the first three components of vector a.

% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

plot3(a(1,:),a(2,:),a(3,:))
hold on
for i = 1:nargin-1
    a = varargin{i};
    plot3(a(1,:),a(2,:),a(3,:))
end
hold off
