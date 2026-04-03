function pp = validateGeo(pp)
% ValidateGeo propagates the geodetic coordinates in the maneuvered and in
% the unmaneuvered case. This allows one to verify if the maneuver is 
% suitable for the respectance of the SK box.
% 
% INPUT:  pp = [struct] Postprocess structure
% 
% OUTPUT: pp = [struct] Postprocess structure
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

if strcmpi(pp.orbit,'geo') && pp.enableSkTarget
    pr              = pp;
    pr.x_s          = pp.validationAbsTraj(:,end);                         % [-] (6,1) Initial state
    pr.dt           = 5000/pp.Tsc;                                         % [-] (1,1) Time step
    pr.T            = pp.T;                                                % [-] (1,1) Orbital period
    dt              = pr.dt;                                               % [-] (1,1)
    pr.tf           = pp.skLength*pr.T;                                             % [-] (1,1) Ending time of simulation
    pr.N            = length(0:dt:pr.tf);                                  % [-] (1,1) Number of nodes
    t_aft           = pp.t(end)+pr.dt:pr.dt:pp.t(end)+(pr.N-1)*pr.dt;      % [-] (1,N) Time discretization
    pp.t_aft        = t_aft;
    pr.et           = pp.et + pp.t(end);                              % [-] (1,1) Ephemeris time
    pr.orbit        = 'geo';
    pr.aida         = pp.aida;
    pr.mu           = pp.mu;
    [~,~,~,lla,~]   = propagateTarget(pr.N,pr.x_s,nan,pr.dt,pr.et,1,pr);   % [-]   (2,N) Propagated geodetic coordinates
    pp.lat_aftMan   = rad2deg(lla(1,2:end));                               % [deg] (1,N) Propagated latitude
    pp.lon_aftMan   = rad2deg(lla(2,2:end));                               % [deg] (1,N) Propagated longitude
    xNoMan          = pp.ballisticTraj(:,end);                             % [-]   (6,1) Unmaneuvered initial state
    [~,~,~,lla,~]   = propagateTarget(pr.N,xNoMan,nan,pr.dt,pr.et,1,pr);   % [-]   (2,N) Unmaneuvered propagated geodetic coordinates
    pp.lat_aftNoMan = rad2deg(lla(1,2:end));                               % [deg] (1,N) Unmaneuvered propagated latitude
    pp.lon_aftNoMan = rad2deg(lla(2,2:end));                               % [deg] (1,N) Unmaneuvered propagated longitude
end
end