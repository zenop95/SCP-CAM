function pp = generateInitShortCanada(pp)
mu     = pp.mu;    % [m^3/s^2]

%% Primary
x0p = toColumn([-1691.2 -3681.3 -5895.7 -4.23 -4.595 4.087]);
primary = cartesian2kepler(x0p,mu);
primary.x0 = x0p;
primary.C0 = zeros(6);
T              = 2*pi/primary.n;         % [s] orbital period
primary.HBR    = 0;           % [km]
primary.mass   = 500;            % [kg] mass
primary.A_drag = 1;              % [m^2] drag surface area
primary.Cd     = 2.2;            % [-] shape coefficient for drag
primary.A_srp  = 1;              % [m^2] SRP surface area
primary.Cr     = 1.31;           % [-] shape coefficient for SRP
primary.cart0 = primary.x0;

%% Secondary
x0s = toColumn([-1691.2+.154 -3681.3 -5895.7 1.074 -2.334 6.524]);
secondary.tca      = 1;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
secondary.x0       = x0s;         
r2e = rtn2eci(x0p(1:3),x0p(4:6));
secondary.relState = [r2e' zeros(3); zeros(3) r2e']*(x0s-x0p);                          % [km] [km/s] Relative cartesian state at TCA
secondary.C0    = eye(6);
secondary.HBR      = 9.3;         % [km]
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