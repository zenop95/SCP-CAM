function pp = generateInitGmm(orbit)

mu  = 398600.4418;    % [km^3/s^2]
C0 = diag(([4.5 10 5 0.15 0.75 0.25]/1e3).^2);
%% Inertial State of the Primary spacecraft
p.C0  = 1/10*C0;             % [km^2] [km^2/s^2] Covariance at TCA
p.inc    = deg2rad(53);            % [rad] inclination
p.a      = 6928; % [km] semimajor axis
p.n      = (mu/p.a^3)^(1/2);    %[rad/s] mean motion
p.ecc      = 0;
p.omega    = 0; 
p.RAAN     = 0;
p.theta0   = 0;            % [rad] inclination
T        = 2*pi/p.n;         % [s] orbital period
p.HBR    = 0.003;           % [km]
p.mass   = 200;            % [kg] mass
p.A_drag = 1;              % [m^2] drag surface area
p.Cd     = 2.2;            % [-] shape coefficient for drag
p.A_srp  = 1;              % [m^2] SRP surface area
p.Cr     = 1.31;           % [-] shape coefficient for SRP
p.cart0  = kepler2cartesian(p.a,p.ecc,p.RAAN,p.inc,p.omega,p.theta0,mu);

%% Relative initial conditions in Hill reference frame
% A0      = 1.000;
% B0      = 0.1;
% y_off   = 0;
% c       = 3;
% alpha0  = deg2rad(15*(c-3));
% beta0   = pi/2 + alpha0;
% n       = p.n;
% r0      = [A0*cos(alpha0); -2*A0*sin(alpha0) + y_off; B0*cos(beta0)];      % [km] Relative position at TCA in LVLH
% v0      = [-n*A0*sin(alpha0); -2*n*A0*cos(alpha0); -n*B0*sin(beta0)];      % [km/s] Relative velocity at TCA in LVLH

%% Relative state of the secondary object at TCA
j             = 1;
s(j).tca      = 1;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
s(j).relState = [0.3 0.3 0.03 1e-4 -1e-4 5e-5]';                        % [km] [km/s] Relative cartesian state at TCA
% s(j).relState = [r0; v0];                        % [km] [km/s] Relative cartesian state at TCA
s(j).x0       = [];
s(j).C0       = 9/10*C0;             % [km^2] [km^2/s^2] Covariance at TCA
s(j).HBR      = p.HBR + 0.035;  % [km]
s(j).mass     = 50;          % [kg] mass
s(j).A_drag   = 1;            % [m^2] drag surface area
s(j).Cd       = 2;          % [-] shape coefficient for drag
s(j).A_srp    = 0.05;            % [m^2] SRP surface area
s(j).Cr       = 1.31;         % [-] shape coefficient for SRP
s(j).x          = [];         % [-] shape coefficient for SRP
s(j).covariance = [];         % [-] shape coefficient for SRP
s(j).w          = 1;
s(j).cdm        = true;
s(j).ang        = false;

pp = struct( ...
            'orbit',      orbit, ...
            'mu',         mu, ...
            'primary',    p, ...
            'secondary',  s, ...
            'T',          T, ....
            'initCovRtn', true ...
            );
end