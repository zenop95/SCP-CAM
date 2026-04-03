function Coe = MEE2COE(Mee)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	 Coe = MEE2COE(Mee)
% 	
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
% 	- a
% 	- e
% 	- in
% 	- OM
% 	- om
% 	- theta
% 	
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p       = Mee(1);
f       = Mee(2);
g       = Mee(3);
h       = Mee(4);
k       = Mee(5);
L       = Mee(6);

a       = p / ( 1 - f^2 - g^2);

e       = sqrt( f^2 + g^2 );

in      = atan2(2 * sqrt(h^2+k^2), 1 - h^2 - k^2);

Om      = atan2(k,h);

om      = atan2(g * h - f * k, f * h + g * k);
if om < 0
    om = 2 * pi + om;
end
    
theta   = L - (Om + om);

Coe = [a; e; in; Om; om; theta];
