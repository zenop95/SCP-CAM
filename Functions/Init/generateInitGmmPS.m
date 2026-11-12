function pp = generateInitGmmPS(pp)

%% Inertial State of the Primary spacecraft
k           = 1;
% p(k).tca    = 1;

p(k).ecc    = 0;              % [-] eccentricity 
p(k).theta0 = 0;              % [rad] initial true anomaly
p(k).omega  = 0;              % [rad] argument of periapsis
p(k).RAAN   = 0;              % [rad] right ascension node longitude
p(k).inc    = 0;              % [rad] inclination
p(k).a      = 6800;           % [km] semimajor axis
p(k).cart0  = kepler2cartesian(p(k).a,p(k).ecc,p(k).RAAN,p(k).inc,p(k).omega,p(k).theta0,pp.mu);

p(k).C0     = diag([0.625 10 3.025 0.00625 0.05625 0.00225])/1e6;  % [km^2] [km^2/s^2] Covariance at TCA
p(k).HBR    = 0;                % [km]
p(k).mass   = 200;              % [kg] mass

p(k).n      = (pp.mu/p(k).a^3)^(1/2);    %[rad/s] mean motion
T           = 2*pi/p(k).n;            % [s] orbital period

p(k).w          = 1;
p(k).cdm        = true;
p(k).ang        = false;

%% Relative state of the secondary object at TCA
j = 1;
s(j).tca      = 1;                                                    % [s] TCA of conjunction w.r.t. initial time t0 = 0
s(j).coe = [6802 4.42e-4 deg2rad(8.4e-5) 0 deg2rad(1.9103) -deg2rad(1.9103)]';
s(j).x0 = kepler2cartesian(s(j).coe(1),s(j).coe(2), s(j).coe(4),s(j).coe(3),s(j).coe(5),s(j).coe(6),pp.mu);                        % [km] [km/s] Relative cartesian state at TCA
s(j).C0       = diag([5.625 90 27.225 0.05625 0.50625 0.02025])/1e6;  % [km^2] [km^2/s^2] Covariance at TCA
s(j).HBR      = p(k).HBR + 0.032;                                     % [km]
s(j).mass     = 50;                                                   % [kg] mass

s(j).w          = 1;
s(j).cdm        = true;
s(j).ang        = false;

pp.primary       = p;
pp.secondary     = s;
pp.T             = T; % Primary's orbital period

% pp.initCovRtn    = true;
end





