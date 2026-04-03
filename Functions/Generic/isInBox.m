function res = isInBox(x0,y0,dx,dy,p)
%isInBox yields 1 if the variable is inside the box defined by the limits
% xlim = [x0 - dx, x0 + dx], ylim = [y0 - dy, y0 + dy]
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
x = p(1);
y = p(2);
res = 1;
if x > x0 + dx || x < x0 - dx || y > y0 + dy || y < y0 - dy
    res = 0;
end

end