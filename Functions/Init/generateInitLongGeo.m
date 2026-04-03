function pp = generateInitLongGeo(pp)

mu  = 398600.4418;    % [km^3/s^2]

%% Inertial state of the Primary
r0 = [9334.75561447905; 41120.9856893676; 1.57083557118839]; % From corrected Laura's case 1
v0 = [-2.99818420674816; .680764971686596; .00637727006252368];
cart = [r0; v0];
p = cartesian2kepler(cart,mu);
p.cart0 = cart;
p.x0    = cart;
p.C0     = diag([0.625 10 3.025 0.00625 0.05625 0.00225])/1e6;
% p.C0     = 0*diag([.45 1 .05 0.0015 0.0075 0.0025]).^2;  % [km^2] [km^2/s^2] Covariance at TCA
p.HBR     = 0.035;  % [km^2] [km^2/s^2] Covariance at TCA
p.mass   = 500;            % [kg] mass
p.A_drag = 1;              % [m^2] drag surface area
p.Cd     = 2.2;            % [-] shape coefficient for drag
p.A_srp  = 1;              % [m^2] SRP surface area
p.Cr     = 1.31;           % [-] shape coefficient for SRP
T = 2*pi/p.n;
%% Relative initial conditions in Hill reference frame
A0      = 1.05;
B0      = 0;
y_off   = 0;
c       = 0;
alpha0  = pi/180*15*(c-3);
beta0   = pi/2 + alpha0;
n       = p.n;
dr0      = [A0*cos(alpha0); -2*A0*sin(alpha0) + y_off; B0*cos(beta0)];     % [km] Relative position at TCA in LVLH
dv0      = [-n*A0*sin(alpha0); -2*n*A0*cos(alpha0); -n*B0*sin(beta0)];     % [km/s] Relative velocity at TCA in LVLH
% dr0      = [1; 0; 0];     % [km] Relative position at TCA in RTN
% dv0      = [0; -1e-3; 0];     % [km/s] Relative velocity at TCA in RTN

%% Relative state of the secondary object at TCA
j = 1;
s(j).tca      = 1;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
s(j).x0       = [];
s(j).relState = [];
s(j).relState = [dr0; dv0];                                                % [km] [km/s] Relative cartesian state at TCA in RTN
% s(j).coe = [42167.76 9e-5 deg2rad(0.119) deg2rad(76.18) deg2rad(241.58) deg2rad(119.44)]';
% s(j).x0 = kepler2cartesian(s(j).coe(1),s(j).coe(2), ...
%     s(j).coe(4),s(j).coe(3),s(j).coe(5),s(j).coe(6),mu);                        % [km] [km/s] Relative cartesian state at TCA
% s(j).C0       = diag([10 5 0.5 0.25 0.75 0.05]/1e3).^2;             % [km^2] [km^2/s^2] Covariance at TCA
s(j).C0       = diag([5.625 90 27.225 0.05625 0.50625 0.02025])/1e6;             % [km^2] [km^2/s^2] Covariance at TCA
s(j).HBR      = p.HBR + 0.010;  % [km]
s(j).mass     = 500;          % [kg] mass
s(j).A_drag   = 1;            % [m^2] drag surface area
s(j).Cd       = 2.2;          % [-] shape coefficient for drag
s(j).A_srp    = 1;            % [m^2] SRP surface area
s(j).Cr       = 1.31;         % [-] shape coefficient for SRP
s(j).x          = [];         % [-] shape coefficient for SRP
s(j).covariance = [];         % [-] shape coefficient for SRP
s(j).w     = 1;
s(j).cdm   = true;
s(j).ang   = false;

pp.fastEncounter = false;                                          
pp.primary       = p;
pp.secondary     = s;
pp.T             = T;
pp.initCovRtn    = true;
pp.singleObject    = true;

end