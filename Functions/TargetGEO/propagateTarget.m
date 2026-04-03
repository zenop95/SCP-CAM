function [dynMaps,llMaps,state,ll_const,nu] = ...
                            propagateTarget(N,x0,newTraj,dt,et,iter,pp)
% InitPropagatorTarget calls the relevant C++ function to perform the
% second order DA propagation that gives the constant part of the
% propagation, the linear maps and the nonlinearity parameter nu.
%     
% INPUT: N            = [-] (1,1) Number of nodes in the optimization
%        newTraj      = [-] (6N,1) Optimization vector resulting from the 
%                       previous iteratoin
%        dt           = [-] (1,1) Time step
%        et           = [-] (1,1) Initial ephemeris time of propagation
%        pp           = [struct] Postprocess structure
% 
% OUTPUT: dynMaps     = [-] (6,6,N) Linear maps for the dynamics
%         llMaps      = [-] (2,3,N) Linear maps from ECI to geodetic   
%         state       = [-] (6,N) constant part of the propagated state
%         ll_const    = [-] (2,N) constant part of the geodetic coordinates
%         nu          = [-] (6,N) nonlinearity parameter nu
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

outDir = "./data_sharing";
if ~exist(outDir,'dir'); mkdir(outDir); end

% ------------------------------------------
% Build JSON input for C++
% ------------------------------------------
input = struct();
input.N       = int64(N);
input.scaling = pp.scaling(1);
input.dt      = dt;
input.et0     = et;

if iter==1
    input.flagProp = int64(1);
    input.x0 = x0(:).';
else
    input.flagProp = int64(2);
    input.x0 = reshape(newTraj,6,[])';   % N×6
end

% Write JSON
fid = fopen(fullfile(outDir,'initial_state.json'),'w');
fwrite(fid,jsonencode(input,'PrettyPrint',true),'char');
fclose(fid);

% ------------------------------------------
% Write AIDA_init.json
% ------------------------------------------
AIDA = struct();
AIDA.flag1   = int64(pp.aida.flag1);
AIDA.flag2   = int64(pp.aida.flag2);
AIDA.flag3   = int64(pp.aida.flag3);
AIDA.gravOrd = int64(pp.aida.gravOrd);
AIDA.mass    = pp.primary.mass;
AIDA.A_drag  = pp.primary.A_drag;
AIDA.Cd      = pp.primary.Cd;
AIDA.A_srp   = pp.primary.A_srp;
AIDA.Cr      = pp.primary.Cr;

fid = fopen(fullfile(outDir,'AIDA_init.json'),'w');
fwrite(fid,jsonencode(AIDA,'PrettyPrint',true),'char');
fclose(fid);

% ------------------------------------------
% Run the C++ executable
% ------------------------------------------
!wsl ./build/bin/findGeo

% ------------------------------------------
% Read output JSON
% ------------------------------------------
output = jsondecode(fileread(fullfile(outDir,'geoOut.json')));

% Output
state    = output.state';
ll_const = output.llConst';
nu       = output.nu';
dynMaps  = permute(output.dynMaps,[2,3,1]);
llMaps   = permute(output.llMaps,[2,3,1]);

end