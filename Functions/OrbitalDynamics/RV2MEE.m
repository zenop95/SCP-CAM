function Mee = RV2MEE(rv, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	rv = MEE2RV(Mee, mu)
% 	
% 	Convert from Modified equinotila elemetns to cartesian coordinates
% 	
% 	INPUT:
% 	- rv 		(1x6) state array (L, L/s)
% 	- mu 		gravitational parameter (L^3/s^2)
%
% 	OUTPUT:
% 	- coe 		(1x6) orbit element array (L, -, rad x 4)
%     .p       
%     .f       
%     .g      
%     .h    
%     .L    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 1
    mu = varargin{1};
else
    mu = 1;
end

Coe = RV2COE(rv, mu);

Mee = COE2MEE(Coe, mu);

