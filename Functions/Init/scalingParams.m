function pp = scalingParams(pp)
% scalingParams
% Computes scaling constants and applies scaling to states, covariance,
% time variables, thrust limits, alt/geo boxes, etc.

% Scaling constants
Lsc = pp.primary.a;                   % [km]     reference length
Vsc = sqrt(pp.mu / Lsc);              % [km/s]   reference velocity
Tsc = Lsc / Vsc;                      % [s]      reference time
Asc = Vsc / Tsc;                      % [km/s^2] reference acceleration

pp.scaling = [Lsc*ones(3,1); Vsc*ones(3,1); Asc*ones(3,1)];
pp.Tsc     = Tsc;

% Covariance transform (ECI/RTN handling if present)
D = diag(1 ./ pp.scaling(1:6));

% Primary objects
for j = 1:numel(pp.primary)
    pp.primary(j).HBR = pp.primary(j).HBR / Lsc;
    pp.primary(j).cart0 = pp.primary(j).cart0 ./ pp.scaling(1:6);
    pp.primary(j).x0    = pp.cart2x(pp.primary(j).cart0);
    
    [r2e, w] = rtn2eci(pp.primary(j).cart0(1:3), pp.primary(j).cart0(4:6));
    R2E      = rot6(r2e, w);
    pp.primary(j).C0 = R2E * D * pp.primary(j).C0 * D' * R2E';           % In ECI Frame
end

% Secondary objects
for j = 1:numel(pp.secondary)
    pp.secondary(j).HBR = pp.secondary(j).HBR / Lsc;
    pp.secondary(j).x0 = pp.secondary(j).x0 ./ pp.scaling(1:6);
    [r2e_s, w_s] = rtn2eci(pp.secondary(j).x0(1:3), pp.secondary(j).x0(4:6));
    R2E_s      = rot6(r2e_s, w_s);
    pp.secondary(j).C0  = R2E_s * D * pp.secondary(j).C0 * D' * R2E_s';  % In ECI Frame
end
% Time variables
pp.dt  = pp.dt_raw / Tsc;
pp.t   = pp.t / Tsc; 
pp.T   = pp.T / Tsc; % Time period of primary

% Epoch (scaled ET at TCA)
etTca  = utc2et(pp.utc);       % requires SPICE helper in path
pp.et  = etTca / pp.Tsc;

% Thrust limits (scaled)
pp.uMax = pp.uMax_unscaled / Asc;
pp.uMin = pp.uMin_unscaled / Asc;

% Objective limit
if strcmpi(pp.obj, 'miss_distance')
    pp.lim = pp.mdLim / Lsc;
else
    pp.lim = pp.pocLim;
end

% Weights finalization
pp.ctrlWeight = 1 / pp.N;  % normalized by number of nodes

% Optional: call switchObj if available
if exist('switchObj','file') == 2
    pp = switchObj(pp);
end

end