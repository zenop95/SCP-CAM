function varargout = ipcErrors(benchmarkIpc,varargin)
% IPCWRAPPER Computes the required Instantaneous Probability of Collision
% metrics, according to the user's input.

% INPUT: numIpc = [-] Precise value of the IPC used as benchmark  
%        Requires at least one additional argument, which needs to be a
%        double scalar

% OUTPUT: variable number of output

% USAGE: [maxErr,cubeErr,constErr] = ipcErr(numIpc,maxIpc,cubeIpc,constIpc);

% Author: Zeno Pavanello, 2022
%--------------------------------------------------------------------------
validateattributes(benchmarkIpc,{'double'},{'vector'})

if nargout == 0 || nargin < 2
    error('ipcError needs at least one IPC input and one output.')
end

benchmarkIpc = toColumn(benchmarkIpc);
N            = length(benchmarkIpc);
varargout    = cell(nargout,1);
for i = 1:nargout
    ipc = varargin{i};
    validateattributes(ipc,{'double'},{'vector','numel',N})
    ipc = toColumn(ipc);
    if N ~= length(ipc)
        error('The two IPC input need to have the same length')
    end
    for j = 1:size(benchmarkIpc)
        aa = benchmarkIpc(benchmarkIpc>0);
        bb = ipc(benchmarkIpc>0);
    end
    absErr = sum(abs(aa-bb))/N;
    relErr = sum(abs(aa-bb)./aa)/N;
    Error = struct('Rel',relErr,'Abs',absErr);
    varargout{i} = Error;
end
    
end

