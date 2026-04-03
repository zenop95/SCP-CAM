function pp = generateInitMultipleStarlink(pp)

mu  = 398600.4418;    % [km^3/s^2]

% Primary spacecraft
primary.C0     = 0*diag([0.5 2 0.1 0.15 0.75 0.25]);  % [km^2] [km^2/s^2] Covariance at TCA
primary.e      = 0;              % [-] eccentricity 
primary.theta0 = 0;              % [rad] initial true anomaly
primary.omega  = 0;              % [rad] argument of periapsis
primary.RAAN   = 0;              % [rad] right ascension node longitude
primary.inc    = deg2rad(53);    % [rad] inclination
primary.a      = 6928;           % [km] semimajor axis
primary.n      = (mu/primary.a^3)^(1/2);    %[rad/s] mean motion
T              = 2*pi/primary.n;         % [s] orbital period
primary.HBR    = 0.003;           % [km]
primary.mass   = 260;            % [kg] mass
primary.A_drag = 1;              % [m^2] drag surface area
primary.Cd     = 2.2;            % [-] shape coefficient for drag
primary.A_srp  = 1;              % [m^2] SRP surface area
primary.Cr     = 1.31;           % [-] shape coefficient for SRP
primary.cart0  = kepler2cartesian(primary.a,primary.e,primary.RAAN,primary.inc,primary.omega,primary.theta0,mu);

%% Secondary structure
secondary = struct();
velAng    = deg2rad([40, 62, 12, 100, 18, 90, 94, 20, 143, 78]); % angle of rotation between primary and secondary velocity. The axis of rotation is the radial direction
angBool   = true;
relState = [0.01 8e-3 -0.02 -0.1 11 0.1; 
            0.01 0.1 0.01 -0.1 0.1 11;
            0 0.05 -0.02 -0  13 0.1; 
            -0.01 0.01 0.07 -0.1 10 0.1; 
            -0.02 0 0.01 2 0 -7;
            0.01 8e-3 -0.02 -0.1 11 0.1; 
            0.01 0 0.01 -0.1 0.1 11;
            0 0.005 -0.02 -0  13 0.1; 
            -0.01 0.01 0.007 -0.1 10 0.1; 
            -0.02 0 0.01 2 0 -7]';
cov      = load('covsStarlink').cov;
tca      = [59; 78; 140; 180; 222; 353; 444; 500; 556; 600];
% tca = tca([1,10]);
% velAng = velAng([1,10]);
% relState = relState(:,[1,10]);
% cov      = cov(:,:,[1,10]);
% cov(:,:,1) = [0.0078      0.0466      0.0644;
%               0.0466      0.2887      0.0103;
%               0.0644      0.0103      0.0145];
% cov(:,:,2) = [ 0.0662  -0.0918  -0.1220;
%               -0.0918   0.1279   0.0879;
%               -0.1220   0.0879   0.1170];
for ii = 1:5
    j = ii;
    secondary(j).tca      = tca(ii);                                      % [s] TCA of conjunction w.r.t. initial time t0 = 0
    secondary(j).x0       = [];         
    secondary(j).mass     = 260;          % [kg] mass
    secondary(j).relState = relState(:,ii);                                 % [km] [km/s] Relative cartesian state at TCA
    secondary(j).velAng   = velAng(ii);                                     % [rad] angle between velocities
    secondary(j).ang      = angBool;                                       % [bool] angle between velocities toggle
    secondary(j).C0       = [cov(:,:,ii) zeros(3,3);zeros(3,6)];            % [km^2] [km^2/s^2] Covariance at TCA
    secondary(j).HBR      = primary.HBR + 0.003;         % [km]
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