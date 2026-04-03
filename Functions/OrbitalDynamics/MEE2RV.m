function rv = MEE2RV(Mee, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	rv = MEE2RV(Mee, mu)
% 	
% 	Convert from Modified equinotila elemetns to cartesian coordinates
% 	
% 	INPUT:
% 	- coe 		(1x6) orbit element array (L, -, rad x 4)
%     .p       
%     .f       
%     .g      
%     .h    
%     .L    
% 	- mu 		gravitational parameter (L^3/s^2)
% 	
% 	OUTPUT:
% 	- rv 		(1x6) state array (L, L/s)
% 	
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 1
    mu = varargin{1};
else
    mu = 1;
end

p       = Mee(1);
f       = Mee(2);
g       = Mee(3);
h       = Mee(4);
k       = Mee(5);
L       = Mee(6);

q           = 1 + f * cos(L) + g * sin(L);
r           = p / q;
alpha2      = h^2 - k^2;
chi2        = h^2 + k^2;
s2          = 1 + chi2;

r_x = r / s2 * ( cos(L) + alpha2 * cos(L) + 2 * h * k * sin(L) );
r_y = r / s2 * ( sin(L) - alpha2 * sin(L) + 2 * h * k * cos(L) );
r_z = 2 * r / s2 * ( h * sin(L) - k * cos(L) );
v_x = - 1 / s2 * sqrt(mu/p) * ( sin(L) + alpha2 * sin(L) - 2 * h * k * cos(L) + g - 2 * f * h * k + alpha2 * g );
v_y = - 1 / s2 * sqrt(mu/p) * ( -cos(L) + alpha2 * cos(L) + 2 * h * k * sin(L) - f + 2 * g * h * k + alpha2 * f );
v_z = 2 / s2 * sqrt(mu/p) * ( h * cos(L) + k * sin(L) + f * h + g * k);

rv = [r_x; r_y; r_z; v_x; v_y; v_z];