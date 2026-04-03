function plot3dbox(x1,x2)
%plot3dbox plots a box given the coordinates of the angles.
% INPUT: 
%        x1 = [3x1] vector of the coordinates of the lowest point
%        x2 = [3x1] vector of the coordinates of the highest point

% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

a = x1(1); b = x1(2); c = x1(3); d = x2(1); e = x2(2); f = x2(3);
plot3([a a], [b b], [c f], 'k--')
hold on
plot3([a a], [b e], [c c], 'k--')
plot3([a d], [b b], [c c], 'k--')
plot3([d d], [e e], [f c], 'k--')
plot3([d d], [e b], [f f], 'k--')
plot3([d a], [e e], [f f], 'k--')
plot3([d d], [e b], [c c], 'k--')
plot3([a d], [e e], [c c], 'k--')
plot3([a a], [b e], [f f], 'k--')
plot3([a a], [e e], [c f], 'k--')
plot3([d d], [b b], [c f], 'k--')
plot3([a d], [b b], [f f], 'k--')
end