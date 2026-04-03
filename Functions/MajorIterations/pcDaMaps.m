function [PoC, pcMaps, xi] = pcDaMaps(relTraj, P, pp)
% pcDaMaps (JSON) computes the linear matrix relating relative position to PoC,
% and returns the constant PoC and trust region, using a JSON-based C++ tool.
%
% INPUT:
%   relTraj : [-] (6, N, M) Relative trajectory of M secondary objects
%   P       : [-] (3, 3, N, M) Covariances at all nodes
%   pp      : [struct] Postprocess structure (expects fields:
%             scaling(1)=Lsc, NCA(1..M), dt, Tsc, obj, PoCType, secondary(i).w,
%             secondary(i).HBR, e2b(3,3,i), timeSubtr, gmmOrder, etc.)
%
% OUTPUT:
%   PoC    : [-] (1,1)  Constant PoC of the mixture
%   pcMaps : [-] (M,3)  Linear maps dPoC/d[r1 r2 r3]_scaled (multiply by Lsc if needed)
%   xi     : [-] (M,3)  Trust region for relative position at conjunction nodes
%
% Notes on scaling:
%   - We pass r (3×1) in scaled units (×Lsc) and P (3×3) in scaled^2 (×Lsc^2)
%   - The C++ returns pcMaps as derivatives w.r.t. the **scaled** coordinates.
%     To get derivatives w.r.t. physical units, multiply by Lsc (we do that here).
%
% JSON files:
%   ./data_sharing/ipcIn.json   (input)
%   ./data_sharing/pocOut.json  (output)
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------   
M   = numel(pp.secondary);
Lsc = pp.scaling(1);

% --------------------------
% Select objective type code
% --------------------------
switch pp.obj
    case 'risk'
        switch lower(pp.PoCType)
            case 'constant', type = 0;
            case 'maximum',  type = 1;
            case 'chan',     type = 2;
            case 'alfano',   type = 3;
            otherwise, error('pcDaMaps:BadType','Invalid PoCType for risk.');
        end
    case 'miss_distance'
        warning('Using linear maps for miss distance is not advisable.'); %#ok<WNTAG>
        type = 4;
    otherwise
        error('pcDaMaps:BadObj','Undefined or invalid PoC model.');
end

% --------------------------
% Build JSON input input
% --------------------------
b1 = tic;
outDir = "./data_sharing";
if ~exist(outDir, 'dir'); mkdir(outDir); end

weights = zeros(M,1);
R       = zeros(M,1);          % HBR * Lsc
r_all   = zeros(M,3);          % r at node NCA(i), scaled by Lsc
P_all   = zeros(M,3,3);        % 3x3 scaled covariances
e2b_all = zeros(M,3,3);        % 3x3 transforms

for i = 1:M
    weights(i) = pp.secondary(i).w;
    R(i)       = pp.secondary(i).HBR * Lsc;

    k = pp.NCA(i);  % conjunction node for i-th secondary
    % relTraj is expected as (6, N, M); use components 1..3
    r_all(i,:) = (reshape(relTraj(1:3, k, i),1,3)) * Lsc;

    % covariance (3x3) scaled^2
    cov_ik = P(:,:,k,i) * (Lsc^2);
    P_all(i,:,:) = cov_ik;

    % direction cosine matrix e2b (3x3)
    e2b_all(i,:,:) = pp.e2b(:,:,i);
end

% P and e2b as list-of-matrices (even when M=1)
P_cells = cell(1,M);
e2b_cells = cell(1,M);

for i = 1:M
    P_cells{i}   = squeeze(P_all(i,:,:));
    e2b_cells{i} = squeeze(e2b_all(i,:,:));
end

input = struct( ...
    'M',        int64(M), ...
    'type',     int64(type), ...
    'weights',  weights(:).', ...
    'R',        R(:).', ...
    'r',        r_all, ...
    'P',        P_cells, ...
    'e2b',      e2b_cells ...
);


% jsonencode supports nested numeric arrays; ensure pretty for debug
fid = fopen(fullfile(outDir, 'ipcIn.json'), 'w');
if fid < 0
    error('pcDaMaps:IO','Cannot open ipcIn.json for writing.');
end
fwrite(fid, jsonencode(input, 'PrettyPrint', true), 'char');
fclose(fid);
a1 = toc(b1);

% --------------------------
% Run C++ executable (JSON)
% --------------------------
!wsl ./build/bin/pocMaps

% --------------------------
% Read JSON output
% --------------------------
b = tic;
raw = fileread(fullfile(outDir, 'pocOut.json'));
res = jsondecode(raw);

% PoC constant
PoC = res.PoC;

% pcMaps: returned as Mx3 derivatives wrt scaled coords -> multiply by Lsc
pcMaps = res.pcMaps * Lsc;

% trustRegion: returned as flat 3*M vector -> reshape to Mx3 and scale
tr = res.trustRegion(:);
if numel(tr) ~= 3*M
    error('pcDaMaps:BadTrustRegion','Expected trustRegion length 3*M.');
end
xi = reshape(tr, 3, M).';   % (M,3)
xi = xi * Lsc;
xi(xi < 1e-2) = 1e-2;

% Accumulate timing: C++ time + MATLAB write/read overhead
pp.timeSubtr = pp.timeSubtr + res.timeMs/1000 + toc(b) + a1;
end