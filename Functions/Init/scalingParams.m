function pp = scalingParams(pp)
% scalingParams
% Computes scaling constants and applies scaling to states, covariance,
% time variables, thrust limits, alt/geo boxes, etc.

% --- Scaling constants
Lsc = pp.primary.a;                   % [km]     reference length
Vsc = sqrt(pp.mu / Lsc);              % [km/s]   reference velocity
Tsc = Lsc / Vsc;                      % [s]      reference time
Asc = Vsc / Tsc;                      % [km/s^2] reference acceleration

pp.scaling = [Lsc*ones(3,1); Vsc*ones(3,1); Asc*ones(3,1)];
pp.Tsc     = Tsc;

% --- Primary
pp.primary.cart0 = pp.primary.cart0 ./ pp.scaling(1:6);
pp.primary.x0    = pp.cart2x(pp.primary.cart0);

% Covariance transform (ECI/RTN handling if present)
D = diag(1 ./ pp.scaling(1:6));
if isfield(pp, 'initCovRtn') && pp.initCovRtn
    [r2e, w] = rtn2eci(pp.primary.cart0(1:3), pp.primary.cart0(4:6));
    R2E      = rot6(r2e, w);
else
    R2E = eye(6);
end
pp.primary.C0 = R2E * D * pp.primary.C0 * D * R2E';

% --- Secondary objects
for j = 1:numel(pp.secondary)
    pp.secondary(j).HBR = pp.secondary(j).HBR / Lsc;
    pp.secondary(j).C0  = D * pp.secondary(j).C0 * D;

    if isfield(pp.secondary(j),'relState') && ~isempty(pp.secondary(j).relState)
        pp.secondary(j).relState = pp.secondary(j).relState ./ pp.scaling(1:6);
    end
    if isfield(pp.secondary(j),'x0') && ~isempty(pp.secondary(j).x0)
        pp.secondary(j).x0 = pp.secondary(j).x0 ./ pp.scaling(1:6);
    end
end

% --- Time variables
pp.dt  = pp.dt_raw / Tsc;
pp.t   = pp.t / Tsc;
pp.T   = pp.T / Tsc;

% Epoch (scaled ET at TCA)
etTca  = utc2et(pp.utc);       % requires SPICE helper in path
pp.et  = etTca / pp.Tsc;

% --- Thrust limits (scaled)
pp.uMax = pp.uMax_unscaled / Asc;
pp.uMin = pp.uMin_unscaled / Asc;

% --- Objective limit
if strcmpi(pp.obj, 'miss_distance')
    pp.lim = pp.mdLim / Lsc;
else
    pp.lim = pp.pocLim;
end

% --- Station-keeping geometry (scaled, rad)
pp.altLim = pp.altLimKm(:) / Lsc;        % [ - ]

% --- Weights finalization
pp.ctrlWeight = 1 / pp.N;                % normalized by number of nodes

% --- Optional: call switchObj if available
if exist('switchObj','file') == 2
    pp = switchObj(pp);
end

end