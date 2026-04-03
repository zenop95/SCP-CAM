function varargout = ipcWrapper(mu,P,R,varargin)
% IPCWRAPPER Computes the required Instantaneous Probability of Collision
% metrics, according to the user's input.

% INPUT: mu = [m] expected value of the relative position random variable
%        P  = [m^2] Covariance Matrix of the relative position random variable
%        R  = [m] Hard Body Radius (HBR) of the collision
%        Requires at least one additional argument, which needs to be a
%        string belonging to the following values
%        accepted values = {'max','cuboid','spheroid','constant','precise',...
%                           'xy','xz','yz'}

% OUTPUT: variable number of output

% USAGE: [maxIpc,cubeIpc,constIpc,numIpc] = ...
%             ipcWrapper(mu,P,R, 'max','cuboid','spheroid','constant',...
%                       'precise','xy','xz','yz');

% Author: Zeno Pavanello, 2022
%--------------------------------------------------------------------------

if nargout == 0 || nargin < 4
    error('ipcWrapper needs at least one IPC type input and one output.')
end

varargout = cell(nargout,1);
for i = 1:nargout
    if strcmpi("max",varargin{i})
        ipc = maximumIpc(mu,P,R);
    elseif strcmpi("cuboid",varargin{i})
        ipc = cuboidIpc(mu,P,R); 
    elseif strcmpi("spheroid",varargin{i})
        ipc = spheroidIpc(mu,P,R); 
    elseif strcmpi("constant",varargin{i})
        ipc = constantIpc(mu,P,R);
    elseif strcmpi("precise",varargin{i})
        ipc = numericIpc(mu,P,R);
    elseif strcmpi("cube precise",varargin{i})
        ipc = cubeNumericIpc(mu,P,R);
    elseif strcmpi("xy",varargin{i})
        ipc = projectedIpc(mu,P,R,'xy');
    elseif strcmpi("xz",varargin{i})
        ipc = projectedIpc(mu,P,R,'xz');
    elseif strcmpi("yz",varargin{i})
        ipc = projectedIpc(mu,P,R,'yz');
    else
        error('Invalid IPC type')
    end
        varargout{i} = ipc;
end

end

