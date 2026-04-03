function pp = generateInitMultiple(orbit)

mu  = 398600.4418;    % [km^3/s^2]

% Primary spacecraft
primary.C0     = 0*diag([0.5 2 0.1 0.15 0.75 0.25]);  % [km^2] [km^2/s^2] Covariance at TCA
primary.e      = 0;              % [-] eccentricity 
primary.theta0 = 0;              % [rad] initial true anomaly
primary.omega  = 0;              % [rad] argument of periapsis
primary.RAAN   = 0;              % [rad] right ascension node longitude
primary.inc    = 1.6;            % [rad] inclination
primary.n      = 1.125915e-3;    %[rad/s] mean motion
primary.a      = (mu/primary.n^2)^(1/3); % [km] semimajor axis
T              = 2*pi/primary.n;         % [s] orbital period
primary.HBR    = 0.02;           % [km]
primary.mass   = 500;            % [kg] mass
primary.A_drag = 1;              % [m^2] drag surface area
primary.Cd     = 2.2;            % [-] shape coefficient for drag
primary.A_srp  = 1;              % [m^2] SRP surface area
primary.Cr     = 1.31;           % [-] shape coefficient for SRP
primary.cart0  = kepler2cartesian(primary.a,primary.e,primary.RAAN,primary.inc,primary.omega,primary.theta0,mu);

%% Secondary structure
secondary = struct();
relState = [0.01 8e-3 -0.02 -0.1 11 0.1; 
            0.01 0 0.01 -0.1 0.1 11;
            0.0 0.005 -0.02 -0  13 0.1; 
            -0.01 0.01 0.007 -0.1 10 0.1; 
            -0.02 0 0.01 2 0 -7]';
cov      = [4.5 10 5 0.15 0.75 0.25;
            2.5 12 10 0.85 0 0.25;
            6 6 5 0.15 0.75 0.25;
            4.5 7 2 0.15 0.75 0.25;
            4.5 10 5 0.15 0.75 0.25]';
tca      = [38; 78; 180; 222; 300];
for j = 1:2
    secondary(j).tca      = tca(j);                                      % [s] TCA of conjunction w.r.t. initial time t0 = 0
    secondary(j).x0       = [];         
    secondary(j).mass     = 100;          % [kg] mass
    secondary(j).relState = relState(:,j);                                 % [km] [km/s] Relative cartesian state at TCA
    secondary(j).C0       = diag((cov(:,j)/1e3).^2);      % [km^2] [km^2/s^2] Covariance at TCA
    secondary(j).HBR      = primary.HBR + 0.01;         % [km]
    secondary(j).A_drag   = 1;            % [m^2] drag surface area
    secondary(j).Cd       = 2.2;          % [-] shape coefficient for drag
    secondary(j).A_srp    = 1;            % [m^2] SRP surface area
    secondary(j).Cr       = 1.31;         % [-] shape coefficient for SRP
    secondary(j).x          = [];         % [-] shape coefficient for SRP
    secondary(j).covariance = [];         % [-] 
    secondary(j).w = 1;        
    secondary(j).cdm = true;        
end
% save('sec','secondary')
% secondary = load('sec').secondary;  
%%

pp.fastEncounter = true;                                          
pp.primary       = primary;
pp.secondary     = secondary;
pp.T             = T;
pp.initCovRtn    = true;
pp.singleObject  = false;

end