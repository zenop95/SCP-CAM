function pp = refineTca(M, t, cart, pp)
% refineTca  Refines the discretized time vector by inserting the real TCAs
%            computed by a JSON-based C++ tool.
%
% INPUTS:
%   M    [-] (1,1)  Number of secondary objects (n_conj * gmmOrder)
%   t    [-] (N,1)  Old time vector
%   cart [-] (6,N)  Cartesian state of the primary at all nodes
%   pp   [struct]   Postprocess structure (expects fields: NCA, dt, Tsc,
%                   gmmOrder, secondary(j).cart(6,N), timeSubtr, etc.)
%
% OUTPUT:
%   pp   [struct]   Updated postprocess structure with refined time grid
%
% Notes:
%   - Writes:  ./data_sharing/initial_state.json
%   - Calls:   !wsl ./build/bin/findRealTca
%   - Reads:   ./data_sharing/tcaOut.json (fields: tca[], timeMs)
%
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

outDir = "./data_sharing";

% Primary / secondary states at the current candidate TCA nodes
% xdump(j,:) = primary state at node NCA(j)
% xdums(j,:) = secondary j state at node NCA(j)
xdump = zeros(6,M);
xdums = zeros(6,M);

for j = 1:M
    k  = pp.NCA(j);
    xdump(:,j) = cart(:, k);                       % primary at node k
    xdums(:,j) = pp.secondary(j).cart(:, k);       % secondary j at k
    relState = xdump(:,j) - xdums(:,j);
    vCond = relState(1:3)'*relState(4:6);
end

input = struct( ...
    'M',     int64(M), ...
    'xdump', xdump', ...
    'xdums', xdums' ...
);

% Write JSON input
fid = fopen(fullfile(outDir, 'initial_state.json'), 'w');
if fid < 0
    error('refineTca:IO', 'Cannot open initial_state.json for writing.');
end
fwrite(fid, jsonencode(input, 'PrettyPrint', true), 'char');
fclose(fid);

% -------------------------
% Run the JSON-enabled TCA finder
% -------------------------
% Assumes your JSON-based binary still uses the same name
% and reads ./data_sharing/initial_state.json
!wsl ./build/bin/findRealTca

% -------------------------
% Read JSON output
% -------------------------
raw = fileread(fullfile(outDir, 'tcaOut.json'));
res = jsondecode(raw);
tcaOff = res.tca(:);     % column vector of size M
timeMs = res.timeMs;     % total timing in milliseconds

% Combine offset with discrete grid indices -> actual TCA absolute time
% (same formula as legacy code)
tca = (pp.NCA(:) - 1 - pp.N_back) * pp.dt + tcaOff;

xpN = propKepOde(xdump(:,4),zeros(3,1),tcaOff(4),1);
xsN = propKepOde(xdums(:,4),zeros(3,1),tcaOff(4),1);
xRelN = xpN-xsN;
res   = xRelN(1:3)'*xRelN(4:6);

% -------------------------
% Introduce real TCAs into the time vector
% -------------------------
n_conj = M / pp.gmmOrder;

% For CDM-constrained conjunctions, remove TCAs (set to NaN)
for j = 1:n_conj
    if pp.secondary(pp.gmmOrder*j).cdm
        idx = (pp.gmmOrder*(j-1)+1) : (pp.gmmOrder*j);
        tca(idx) = nan;
    else
        % Remove TCAs too close to the "central" mixand TCA (< 1 sec)
        cidx = ceil(pp.gmmOrder*( (j-1) + 0.5));   % central index in group
        refTca = tca(cidx);
        % All but the central element
        left  = 1:floor(pp.gmmOrder/2);
        right = ceil(pp.gmmOrder/2)+1 : pp.gmmOrder;
        for jj = [left, right]
            kk = pp.gmmOrder*(j-1) + jj;
            if abs(tca(kk) - refTca) * pp.Tsc < 1
                tca(kk) = nan;
            end
        end
    end
end

% Remove original discrete instants too close to true TCA (< 1 sec)
for j = 1:M
    k = pp.NCA(j);
    if abs(t(k) - tca(j)) * pp.Tsc < 1
        t(k) = nan;
    end
end

% Merge, sort, and drop NaNs
t = sort([t; tca], 'ascend');
t(isnan(t)) = [];

% -------------------------
% Fix simulation parameters to the new time vector
% -------------------------
N = numel(t);

for j = 1:n_conj
    grp = (j-1)*pp.gmmOrder + (1:pp.gmmOrder);

    if ~pp.secondary(grp(end)).cdm
        % Snap each non-NaN TCA to nearest time node index
        for k = grp
            if ~isnan(tca(k))
                [~, ind] = min(abs(t - tca(k)));
                pp.NCA(k) = ind;
            end
        end

        % For NaN TCAs, inherit the central mixand index
        cidx = ceil(pp.gmmOrder*( (j-1) + 0.5 ));
        for k = grp
            if isnan(tca(k))
                pp.NCA(k) = pp.NCA(cidx);
            end
        end
    end
end

pp.N_forw = nnz(t >= 0) - 1;
pp.N_back = N - pp.N_forw - 1;
pp.t      = t(:);
pp.N      = N;

for j = 1:M
    pp.secondary(j).tca = pp.NCA(j) - pp.N_back;
end

% Accumulate timing from the C++ JSON tool
pp.timeSubtr = pp.timeSubtr + timeMs / 1000;

end