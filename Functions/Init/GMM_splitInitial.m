function pp = GMM_splitInitial(pp)
% GMM_splitInitial Splits each primary & secondary into a Gaussian mixture in ECI/RTN frames 
% GMM_splitInitial takes the postprocess structure containing the objects and replaces each primary & secondary 
% by gmmOrder components following a nonlinearity-guided split. For primaries or secondaries given in RTN (relState)
% the function converts to ECI, performs the split, then converts states/covariances back to the requested frame(s).
%
% BEHAVIOR:
%  - For each entry in pp.secondary:
%    * If a secondary provides x0 it is used as ECI state.
%    * Covariances are transformed between RTN and ECI.
%    * A direction of maximum nonlinearity is found by calling an external
%      executable (writes initial_state.json reads trustRegion.json).
%    * The Vittaldev_Multipledir split is applied to produce xi (states), Ci (covariances)
%      and wi (weights) for each Gaussian component; pp.secondary is updated.
%
% INPUT:
%   pp = [struct] Postprocess structure with (used fields)
%        .primary     (1,M) struct array with fields x0, C0, cdm, ...
%        .secondary   (1,M) struct array with fields x0, relState, C0, tca, cdm, ...
%        .gmmOrder    (1,2) number of Gaussian components per primary and secondary
%        .primary.x0 (6,1) primary absolute state used for frame conversion
%        .N, .dt, .et  scalars used to generate initial_state.dat for nliGmm
%
% OUTPUT:
%   pp = [struct] same input struct with pp.primary and pp.secondary expanded/updated:
%        - each original primary replaced by gmmOrder(1) entries
%        - each original secondary replaced by gmmOrder(2) entries
%        - for each component: .w, .x0 (ECI), and .C0 updated
%
% EXAMPLE USAGE:
%   pp = GMM_splitInitial(pp);
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

p = pp.primary;
n1 = pp.gmmOrder(1);

s = pp.secondary;
n2 = pp.gmmOrder(2);

%% Transform primary into ECI
original = p(1);
if ~isempty(p(1).cart0)  
    state = p(1).cart0;   % [-] (6,1) Absolute primary state in ECI at the time of conjunction        
    stateNli  = state;
else
    error('absolute state of the primary must be defined');
end

%% Find direction of maximum nonlinearity
input      = struct();
input.dt   = pp.N*pp.dt; 
input.x0   = stateNli; 

fid = fopen('./data_sharing/nli.json','w'); 
fwrite(fid,jsonencode(input),'char'); 
fclose(fid);

!wsl ./build/bin/nliGmm
nli = jsondecode(fileread('./data_sharing/trustRegion.json')).trustRegion;

%% Perform the split along the found direction
dirRtn = normalize(nli.*diag(C0rtn),'norm');
dir    = R2Ep*dirRtn;
xc     = state;
C0c    = C0;
[xi, Ci, wi] = VittaldevSplit(xc,C0c,dir,n1);
for j = 1:n1
    priInd                  = j + (k-1)*n1;           
    pp.primary(priInd)    = original;
    pp.primary(priInd).w  = wi(j);
    pp.primary(priInd).x0 = xi(:,j);                     % state in ECI
    % if isempty(s(k).x0) && ~isempty(s(k).relState)
    %     pp.secondary(secInd).relState = R2Ep\(pp.secondary(secInd).x0 - pp.primary.cart0);
    %     pp.secondary(secInd).x0 = [];                     % state in ECI
    % end
    if pp.initCovRtn_P
        pp.primary(priInd).C0 = E2Rp*Ci(:,:,j)*E2Rp';        % covariance in RTN of primary
    else
        pp.primary(priInd).C0 = Ci(:,:,j);
    end
end

%% Transform secondary into ECI  
original = s(1);
if ~isempty(s(1).x0) 
    state = s(1).x0;     % [-] (6,1) Absolute secondary state in ECI at the time of conjunction        
    stateNli  = state;
else
    error('The absolute state of the secondary must be defined');
end

%% Find direction of maximum nonlinearity
input      = struct();
input.dt   = pp.N*pp.dt; 
input.x0   = stateNli; 

fid = fopen('./data_sharing/nli.json','w'); 
fwrite(fid,jsonencode(input),'char'); 
fclose(fid);

!wsl ./build/bin/nliGmm
nli = jsondecode(fileread('./data_sharing/trustRegion.json')).trustRegion;

%% Perform the split along the found direction
dirRtn = normalize(nli.*diag(C0rtn),'norm');
dir    = R2Es*dirRtn;
xc     = state;
C0c    = C0;
[xi, Ci, wi] = VittaldevSplit_Multipledir(xc,C0c,dir,n2);
for j = 1:n2
    secInd                  = j + (k-1)*n2;           
    pp.secondary(secInd)    = original;
    pp.secondary(secInd).w  = wi(j);
    pp.secondary(secInd).x0 = xi(:,j);                     % state in ECI
    if isempty(s(k).x0) && ~isempty(s(k).relState)
        pp.secondary(secInd).relState = R2Ep\(pp.secondary(secInd).x0 - pp.primaryReferenceState);
        pp.secondary(secInd).x0 = [];                     % state in ECI
    end
    if pp.initCovRtn_S
        pp.secondary(secInd).C0 = E2Rs*Ci(:,:,j)*E2Rs';        % covariance in RTN of secondary
    else
        pp.secondary(secInd).C0 = Ci(:,:,j);
    end
end
end
