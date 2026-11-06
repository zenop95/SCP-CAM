function pp = generateInitGmmPS(pp)

mu  = 398600.4418;    % [km^3/s^2]
%% Inertial State of the Primary spacecraft
k           = 1;
    p(k).tca    = 1;
p(k).C0     = diag([0.625 10 3.025 0.00625 0.05625 0.00225])/1e6;             % [km^2] [km^2/s^2] Covariance at TCA
p(k).ecc    = 0;              % [-] eccentricity 
p(k).theta0 = 0;              % [rad] initial true anomaly
p(k).omega  = 0;              % [rad] argument of periapsis
p(k).RAAN   = 0;              % [rad] right ascension node longitude
p(k).inc    = 0;    % [rad] inclination
p(k).a      = 6800;           % [km] semimajor axis
p(k).n      = (mu/p(k).a^3)^(1/2);    %[rad/s] mean motion
T           = 2*pi/p(k).n;         % [s] orbital period
% p(k).HBR    = 0.003;          % [km]
p(k).HBR    = 0;                % [km]
p(k).mass   = 200;              % [kg] mass

% Neglect (as AIDA set to 0) 
p(k).A_drag = 1;              % [m^2] drag surface area
p(k).Cd     = 2.2;            % [-] shape coefficient for drag
p(k).A_srp  = 1;              % [m^2] SRP surface area
p(k).Cr     = 1.31;           % [-] shape coefficient for SRP

p(k).cart0  = kepler2cartesian(p(k).a,p(k).ecc,p(k).RAAN,p(k).inc,p(k).omega,p(k).theta0,mu);
   % p(k).x0     = p(k).cart0;
    
    p(k).w          = 1;
    p(k).cdm        = true;
    p(k).ang        = false;

%% Relative initial conditions in Hill reference frame
% A0      = 1.000;
% B0      = 0.1;
% y_off   = 0;
% c       = 5;
% alpha0  = deg2rad(15*(c-3));
% beta0   = pi/2 + alpha0;
% n       = p(k).n;
% r0      = [A0*cos(alpha0); -2*A0*sin(alpha0) + y_off; B0*cos(beta0)];      % [km] Relative position at TCA in LVLH
% v0      = [-n*A0*sin(alpha0); -2*n*A0*cos(alpha0); -n*B0*sin(beta0)];      % [km/s] Relative velocity at TCA in LVLH

%% Relative state of the secondary object at TCA
j = 1;
s(j).tca      = 1;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
       % s(j).relState = [r0; v0];                                                  % [km] [km/s] Relative cartesian state at TCA in RTN
       % s(j).relState = [1 0.1 0 0 -2e-3 -1e-5]';                        % [km] [km/s] Relative cartesian state at TCA
       % s(j).x0       = [];

% Relative state of secondary as expressed in primary's RTN at TCA:
% s(j).relState = [-1.00481459101518
%                          0
%                          0
%      -0.000112790199313706
%        0.00338808867193308
%       1.12279140452717e-05];                                                  % [km] [km/s] Relative cartesian state at TCA in RTN
s(j).relState       = [];
s(j).coe = [6802 4.42e-4 deg2rad(8.4e-5) 0 deg2rad(1.9103) -deg2rad(1.9103)]';
s(j).x0 = kepler2cartesian(s(j).coe(1),s(j).coe(2), ...
    s(j).coe(4),s(j).coe(3),s(j).coe(5),s(j).coe(6),mu);                        % [km] [km/s] Relative cartesian state at TCA
s(j).C0       = diag([5.625 90 27.225 0.05625 0.50625 0.02025])/1e6;             % [km^2] [km^2/s^2] Covariance at TCA
s(j).HBR      = p(k).HBR + 0.032;  % [km]
s(j).mass     = 50;          % [kg] mass

% Neglect (as AIDA set to 0) 
s(j).A_drag   = 1;            % [m^2] drag surface area
s(j).Cd       = 2;          % [-] shape coefficient for drag
s(j).A_srp    = 0.05;            % [m^2] SRP surface area
s(j).Cr       = 1.31;         % [-] shape coefficient for SRP
s(j).x          = [];         % [-] shape coefficient for SRP
s(j).covariance = [];         % [-] shape coefficient for SRP

s(j).w          = 1;
s(j).cdm        = true;
s(j).ang        = false;


pp.fastEncounter = false; 

  pp.primaryReference = p(1);
  pp.primaryReferenceState = p(1).cart0(:);

pp.primary       = p;
pp.secondary     = s;
pp.T             = T;

%pp.initCovRtn    = true;
pp.initCovRtn_P    = true;
pp.initCovRtn_S    = true;

pp.singleObject  = true;
%pp.singleObject_P  = true; % Primary Object
%pp.singleObject_S  = true; % Secondary Object
end





