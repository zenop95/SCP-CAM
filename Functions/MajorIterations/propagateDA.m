function [dynMaps,STM,llMaps,cartMaps,state,cart,lla_const,nu,pp] = ...
                                  propagateDA(iter,oldTraj,validation,pp)
% PropagateDA calls the relevant C++ function to perform the first or
% second order DA propagation that gives the constant part of the
% propagation, the linear maps and the nonlinearity parameter nu.
% 
% INPUT:  t          = [-] (N,1) Time vector
%         iter       = [-] (1,1) Current iteration counter
%         oldTraj    = [-] (6N,1) Relative trajectory from previous iteration
%         validation = [bool] (1,1) Flag to perform validation propagation
%         pp         = postprocess structure
%
% OUTPUT: dynMaps  = [-] (6,9,N) Linear maps for the dynamics
%         STM      = [-] (6,6,N) State transition matrices
%         llMaps   = [-] (2,3,N) Linear maps from state to geodetic   
%         cartMaps = [-] (6,6,N) Linear maps from state to cartesian   
%         state    = [-] (6,N)   Constant part of the propagated state
%         cart     = [-] (6,N)   Constant part of the cartesian coordinates
%         ll_const = [-] (2,N)   Constant part of the geodetic coordinates
%         nu       = [-] (6,N)   Nonlinearity parameter nu
%         pp       = [struct]    Postprocess structure. 
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
N=pp.N; t=pp.t;
if validation 
    order = 1; 
    u     = reshape(pp.u,3*N,1); 
else 
    order = 2; 
    u     = zeros(3*N,1); 
end

input         = struct();
input.N       = N; 
input.t       = t; 
input.et      = pp.et;
input.order   = order; 
input.uMax    = pp.uMax;
input.scaling = pp.scaling(1);
input.flagRtn = pp.flagRtn;

if iter==1||validation
    input.flagProp = 1;
    input.x0       = pp.primary.x0(:)';
    input.u        = reshape(u,3,[])';
else
    input.flagProp = 2;
    input.oldTraj  = reshape(oldTraj,6,[])';
    input.u        = reshape(pp.majorIter(end).u,3,[])';
end

fid = fopen('./data_sharing/initial_state.json','w'); 
fwrite(fid,jsonencode(input),'char'); 
fclose(fid);
!wsl ./build/bin/propAida
output = jsondecode(fileread('./data_sharing/out_prop.json'));

state     = output.state';
cart      = state;
lla_const = output.latlonConst';

% build dynMaps (9x6xN), cartMaps(6x6), llMaps(9x3)
dynMaps  = permute(output.dynMaps,[2,3,1]);          % 9x6xN
llMaps   = permute(output.latlonMaps,[2,3,1]);       % 9x3xN
STM      = dynMaps(1:6,1:6,:);
llMaps   = llMaps(:,1:6,:);
cartMaps = repmat(eye(6),[1,1,N]);
if isfield(output,'trustRegion')
    nu = output.trustRegion'; 
else 
    nu = []; 
end

pp.timeSubtr = pp.timeSubtr+output.timeMs/1000;
end
