function [s, c, r, t] = atan2sinCos(y, x)
% ATAN2SC: sin and cosine of atan(y/x) 
%
% Usage:
%   [s, c, r, t] = atan2sinCos(y, x)
%
% Inputs:
%   y: 1xN fraction numerator
%   x: 1xN fraction denumarator
%
% Outputs:
%    s: 1xN sin(t) where tan(t) = y/x
%    c: 1xN cos(t) where tan(t) = y/x
%    r: 1xN sqrt(x^2 + y^2)
%    t: 1xN arctan of y/x

reshape(y, 1, []);
reshape(x, 1, []);

if ~isequal(size(x, 2), size(y, 2))
   error('''x'' and ''y'' sizes do not match.'); 
end

r = sqrt(y.^2 + x.^2);
s = y./r;
c = x./r;
t = atan2(s, c);
end