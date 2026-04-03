function pp = generateInitShort(pp)
ind    = pp.ESACase;
mu     = pp.mu;    % [m^3/s^2]
load("data/conjunctions_leo.mat");
% load("train.mat");
B = table2array(data); clear data;

%% Primary
x0p = toColumn(B(ind,3:8));
primary = cartesian2kepler(x0p,mu);
primary.x0 = x0p;
primary.C0 = [[B(ind,9) B(ind,12)  B(ind,13);
           B(ind,12)  B(ind, 10) B(ind,14);
           B(ind,13)  B(ind, 14)  B(ind,11)] zeros(3,3); 
                                            zeros(3,6)];
T              = 2*pi/primary.n;         % [s] orbital period
primary.HBR    = B(ind,2)/2;           % [km]
primary.mass   = 500;            % [kg] mass
primary.A_drag = 1;              % [m^2] drag surface area
primary.Cd     = 2.2;            % [-] shape coefficient for drag
primary.A_srp  = 1;              % [m^2] SRP surface area
primary.Cr     = 1.31;           % [-] shape coefficient for SRP
primary.cart0 = primary.x0;

%% Secondary
x0s                = toColumn(B(ind,15:20)); % [km] [km/s] Secondary initial state in ECI
secondary.tca      = 1;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
secondary.x0       = x0s;         
r2e = rtn2eci(x0p(1:3),x0p(4:6));
secondary.relState = [r2e' zeros(3); zeros(3) r2e']*(x0s-x0p);                          % [km] [km/s] Relative cartesian state at TCA
secondary.C0    = [[B(ind,21) B(ind,24) B(ind,25);
           B(ind,24)  B(ind,22) B(ind,26);
           B(ind,25)  B(ind,26) B(ind,23)] zeros(3,3); 
                                            zeros(3,6)];
secondary.HBR      = B(ind,2)/2 + primary.HBR;         % [km]
secondary.mass     = 100;          % [kg] mass
secondary.A_drag   = 1;            % [m^2] drag surface area
secondary.Cd       = 2.2;          % [-] shape coefficient for drag
secondary.A_srp    = 1;            % [m^2] SRP surface area
secondary.Cr       = 1.31;         % [-] shape coefficient for SRP
secondary.x          = [];         % [-] 
secondary.covariance = [];         % [-] 
secondary.w          = 1;        
secondary.cdm        = true;        
secondary.ang        = false;     

pp.fastEncounter = true;                                          
pp.primary       = primary;
pp.secondary     = secondary;
pp.T             = T;
pp.initCovRtn    = true;
pp.singleObject  = true;

end