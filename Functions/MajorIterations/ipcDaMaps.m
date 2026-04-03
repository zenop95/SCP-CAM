function [ipcMaps, xi] = ipcDaMaps(relTraj, P, pp)
% ipcDaMaps (JSON version)
% Computes linear maps relating relative position to IPC and the trust-region xi.
%
% INPUT:
%   relTraj : (6, N, M)
%   P       : (3,3,N,M)
%   pp      : struct with fields:
%             N, secondary(j).w, secondary(j).HBR, scaling(1), NCA, obj, etc.
%
% OUTPUT:
%   ipcMaps : (3, M, N)
%   xi      : (3, M, N)
%
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------   
    

N   = pp.N;
M   = numel(pp.secondary);
Lsc = pp.scaling(1);

% Determine IPC type
switch lower(pp.obj)
    case 'risk'
        type = 0;
    case 'max_risk'
        type = 1;
    case 'd_miss'
        warning("IPC not used for d_miss objective.");
        type = 99;  % not used
    otherwise
        error("Invalid pp.obj");
end

% ---------------------------------------------
% Build JSON input
% ---------------------------------------------
outDir = "./data_sharing";
if ~exist(outDir,'dir'); mkdir(outDir); end

weights = zeros(1,M);
R       = zeros(1,M);

r_all   = zeros(N, M, 3);     % we store coordinates per-node, per-object, 3-dim
P_all   = zeros(N, M, 3, 3);
for j = 1:M
    weights(j) = pp.secondary(j).w;
    R(j)       = pp.secondary(j).HBR * Lsc;
    for i = 1:N
        r_all(i,j,:) = relTraj(1:3,i,j) * Lsc;
        P_all(i,j,:,:) = P(:,:,i,j) * Lsc^2;
    end
end

% Wrap every 3×3 matrix into a cell, to guarantee JSON array
P_cells = cell(N, M);
for i = 1:N
    for j = 1:M
        P_cells{i,j} = squeeze(P_all(i,j,:,:));
    end
end

% JSON payload
payload = struct( ...
    'N',      int64(N), ...
    'M',      int64(M), ...
    'type',   int64(type), ...
    'weights', weights, ...
    'R',       R, ...
    'r',       r_all, ...     % stored as N×M×3
    'P',       P_cells ...    % stored as {N,M} cells → JSON nested arrays
);

% Write JSON
jsonPath = fullfile(outDir, "pocIn.json");
fid = fopen(jsonPath, 'w');
if fid < 0, error("Cannot write ipcIn.json"); end
fwrite(fid, jsonencode(payload, 'PrettyPrint', true), 'char');
fclose(fid);

% ---------------------------------------------
% Run JSON IPC executable
% ---------------------------------------------
!wsl ./build/bin/ipcMaps

% ---------------------------------------------
% Read output JSON
% ---------------------------------------------
raw = fileread(fullfile(outDir,"pocOut.json"));
res = jsondecode(raw);

% Extract IPC maps: size (N × M × 3)
maps = res.ipcMaps;     % as N×M×3 (if JSON produced it that way)

% If decoded as Mx3xN or Nx3xM, normalize:
if ndims(maps) == 3 && size(maps,3)==3 && size(maps,2)==M
    % (N,M,3) -- expected
    ipcMaps = permute(maps, [3 2 1]);   % -> (3,M,N)
elseif ndims(maps)==3 && size(maps,1)==M && size(maps,2)==3 && size(maps,3)==N
    ipcMaps = maps;   % already (M,3,N)
    ipcMaps = permute(ipcMaps, [2 1 3]);   % -> (3,M,N)
else
    error("Unexpected ipcMaps JSON dimensions");
end

ipcMaps = ipcMaps * Lsc;

% Trust region: "trustRegion" is vector of length 3*M*N
tr = res.trustRegion(:);

if numel(tr) ~= 3*M*N
    error("trustRegion length mismatch");
end

xi = reshape(tr, [3, M, N]);
xi = xi * Lsc;
xi(isnan(xi)) = 1e-1;
xi(isinf(xi)) = 100;

% Time accumulation
pp.timeSubtr = pp.timeSubtr + res.timeMs/1000;
end