function pp = simProperties(i)
% simParams
% Assembles the full configuration structure 'pp' by:
%   1) Loading user parameters
%   2) Loading default parameters
%   3) Building the scenario (unscaled)
%   4) Applying scaling and derived limits
%
% Author: Zeno Pavanello
% E-mail: zpav176@aucklanduni.ac.nz
% Date:   2022-2024
%--------------------------------------------------------------------------

% 1. User params
up = userParams(i);

% 2. House defaults
dp = defaultParams();

% 3. Scenario (unscaled)
pp = scenarioParams(up, dp);

% 4. Scaling, limits, and derived values
pp = scalingParams(pp);

% 5. Final overrides and consistency
% If PoC is enforced via constraint -> force minor iter = 1
if pp.pocConstr
    pp.iterMaxMin = 1;
end

% If linear constraint is embedded -> disable pocConstr
pp.embedLinearConstraint = true*pp.pocConstr;
if pp.embedLinearConstraint
    pp.pocConstr = false;
end

% If the gravity order is 0 set it to 1 for incompatibility with AIDA code
if pp.aida.gravOrd == 0
    pp.aida.gravOrd = 1;
end

% If the orbit is GEO, turn off the atmospheric model
if strcmpi(pp.orbit,'geo')
    pp.aida.flag1 = 0;
end

end