function pp = generateInit_CARA(pp,type)
%
% Author: Zeno Pavanello
% E-mail: zpav176@aucklanduni.ac.nz
% Date:   2022-2024
%--------------------------------------------------------------------------
ind    = pp.ESACase;
mu     = pp.mu;    % [m^3/s^2]
load("data/CARA_data.mat");
% load("train.mat");
B = conjunctionData(ind); clear conjunctionData;

%% Primary
x0p            = [B.r1; B.v1]/1e3;
primary        = cartesian2kepler(x0p,mu);
primary.x0     = x0p;
primary.C0     = B.C1/1e6;
T              = 2*pi/primary.n;         % [s] orbital period
primary.HBR    = B.HBR/1e3/2;           % [km]
primary.mass   = 500;            % [kg] mass
primary.A_drag = 1;              % [m^2] drag surface area
primary.Cd     = 2.2;            % [-] shape coefficient for drag
primary.A_srp  = 1;              % [m^2] SRP surface area
primary.Cr     = 1.31;           % [-] shape coefficient for SRP
primary.cart0  = primary.x0;

%% Secondary
x0s                  = [B.r2; B.v2]/1e3; % [km] [km/s] Secondary initial state in ECI
secondary.tca        = 1;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
secondary.x0         = x0s;         
r2e                  = rtn2eci(x0p(1:3),x0p(4:6));
secondary.relState   = [r2e' zeros(3); zeros(3) r2e']*(x0s-x0p);           % [km] [km/s] Relative cartesian state at TCA
secondary.C0         = B.C2/1e6;
secondary.HBR        = B.HBR/1e3/2 + primary.HBR;         % [km]
secondary.mass       = 100;          % [kg] mass
secondary.A_drag     = 1;            % [m^2] drag surface area
secondary.Cd         = 2.2;          % [-] shape coefficient for drag
secondary.A_srp      = 1;            % [m^2] SRP surface area
secondary.Cr         = 1.31;         % [-] shape coefficient for SRP
secondary.x          = [];         % [-] 
secondary.covariance = [];         % [-] 
secondary.w          = 1;        
secondary.cdm        = true;        
secondary.ang        = false;     

if     strcmpi(type,'short'), pp.fastEncounter = true;
elseif strcmpi(type,'long'),  pp.fastEncounter = false; end
pp.primary       = primary;
pp.secondary     = secondary;
pp.T             = T;
pp.initCovRtn    = false;
pp.singleObject  = true;

end