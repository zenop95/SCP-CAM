function Mee = COE2MEE(Coe, mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Mee = COE2MEE(Coe, mu)
% 	
% 	
% 	INPUT:
% 	- Mee 		(1x6) orbit element array (L, -, rad x 4)
%     .p       
%     .f       
%     .g      
%     .h    
%     .L    
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
a       = Coe(1);
e       = Coe(2);
in      = Coe(3);
Om      = Coe(4);
om      = Coe(5);
theta   = Coe(6);

p = a * ( 1 - e^2);
f = e * cos(om + Om);
g = e * sin(om + Om);
h = tan(in/2) * cos(Om);
k = tan(in/2) * sin(Om);
L = Om + om + theta;
L = mod(L,2*pi);
L(L < 0) = L+2*pi;
Mee = [p; f; g; h; k; L];
