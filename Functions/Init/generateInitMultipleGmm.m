function pp = generateInitMultipleGmm(pp)

mu  = 398600.4418;    % [km^3/s^2]

% Primary spacecraft
primary.C0     = 0*diag(([1e-3 1e-2 1e-3 1e-6 5e-6 1e-6]).^2);  % [km^2] [km^2/s^2] Covariance at TCA
primary.e      = 0;              % [-] eccentricity 
primary.theta0 = 0;              % [rad] initial true anomaly
primary.omega  = 0;              % [rad] argument of periapsis
primary.RAAN   = 0;              % [rad] right ascension node longitude
primary.inc    = 1.6;            % [rad] inclination
primary.n      = 1.125915e-3;    %[rad/s] mean motion
primary.a      = (mu/primary.n^2)^(1/3); % [km] semimajor axis
T              = 2*pi/primary.n;         % [s] orbital period
primary.T      = T;             
primary.HBR    = 0.02;           % [km]
primary.mass   = 500;            % [kg] mass
primary.A_drag = 1;              % [m^2] drag surface area
primary.Cd     = 2.2;            % [-] shape coefficient for drag
primary.A_srp  = 1;              % [m^2] SRP surface area
primary.Cr     = 1.31;           % [-] shape coefficient for SRP
primary.cart0  = kepler2cartesian(primary.a,primary.e,primary.RAAN,primary.inc,primary.omega,primary.theta0,mu);

%% Secondary structure
mult = 2;
secondary        = struct();
secondary(1).inc    = -1.6;                                   % [rad] inclination
secondary(1).RAAN   = 0;                                       % [rad] right ascension node longitude
secondary(1).theta0 = 1e-5;                                       % [rad] initial true anomaly
% secondary(1).RAAN   = deg2rad(1.72e-3);                                       % [rad] right ascension node longitude
% secondary(1).theta0 = deg2rad(0.11);                                       % [rad] initial true anomaly
secondary(1).omega  = 0;                                       % [rad] argument of periapsis
secondary(1).T      = mult*primary.T;                          % [s] orbital period
secondary(1).rp     = primary.a;                               % [km] perigee radius
secondary(1).n      = 2*pi/secondary.T;                        %[rad/s] mean motion
secondary(1).a      = (mu/secondary.n^2)^(1/3);                % [km] semimajor axis
secondary(1).e      = 1-secondary.rp/secondary.a;              % [-] eccentricity 
secondary(1).cart0  = kepler2cartesian(secondary.a,secondary.e,secondary.RAAN,secondary.inc,secondary.omega,secondary.theta0,mu);
secondary(1).cdm    = true;

secondary(1).tca      = 1;                                           % [s] TCA of conjunction w.r.t. initial time t0 = 0
secondary(1).x0       = secondary.cart0;         
secondary(1).mass     = 100;          % [kg] mass
secondary(1).relState = [];                        % [km] [km/s] Relative cartesian state at TCA
% secondary(1).C0       = diag([.045 .1 .05 1.5e-5 7.5e-5 2.5e-5]).^2;       % [km^2] [km^2/s^2] Covariance at TCA in secondary RTN
secondary(1).C0       = diag([.075 .1 .0075 1.5e-3 1.5e-3 0.5e-4]).^2;       % [km^2] [km^2/s^2] Covariance at TCA in secondary RTN
secondary(1).HBR      = primary.HBR + 0.01;         % [km]
secondary(1).A_drag   = 1;            % [m^2] drag surface area
secondary(1).Cd       = 2.2;          % [-] shape coefficient for drag
secondary(1).A_srp    = 1;            % [m^2] SRP surface area
secondary(1).Cr       = 1.31;         % [-] shape coefficient for SRP
secondary(1).x          = [];         % [-] shape coefficient for SRP
secondary(1).covariance = [];         % [-] 
secondary(1).w = 1;        

secondary(2)     = secondary(1);
secondary(2).tca = mult*60+1;
secondary(2).cdm = false;

for j = 2:2
    secondary(j) = secondary(1);
    secondary(j).tca = (j-1)*mult*60+1;
    secondary(j).cdm = false;
end

%%


pp.fastEncounter   = true;                                          
pp.primary         = primary;
pp.secondary       = secondary;
pp.T               = T;
pp.initCovRtn      = true;
pp.toggleRefineTca = true;
pp.singleObject    = true;

end