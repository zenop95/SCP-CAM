function pp = splitInitial(pp)
% splitInitial Splits each secondary into a Gaussian mixture in ECI/RTN frames
%
% splitInitial takes the postprocess structure containing secondary objects
% 
% and replaces each secondary by gmmOrder components following a nonlinearity-
% 
% guided split. For secondaries given in RTN (relState) the function converts
% 
% to ECI, performs the split, then converts states/covariances back
% 
% to the requested frame(s).
%
% BEHAVIOR:
%
%  - For each entry in pp.secondary:
% 
%    * If pp.singleObject is true, the function duplicates the first block
%      of gmmOrder components for subsequent objects while copying tca/cdm.
% 
%    * If a secondary provides relState (RTN) it is converted to ECI using
%      pp.primary.cart0 and the primary RTN→ECI DCM.
% 
%    * If a secondary provides x0 it is used as ECI state.
% 
%    * Covariances are transformed between RTN and ECI depending on
%      pp.initCovRtn.
% 
%    * A direction of maximum nonlinearity is found by calling an external
%      executable (writes initial_state.dat and reads trustRegion.dat).
% 
%    * The Vittaldev split is applied to produce xi (states), Ci (covariances)
%      and wi (weights) for each Gaussian component; pp.secondary is updated.
%
% INPUT:
%
%   pp = [struct] Postprocess structure with (used fields)
% 
%        .secondary   (1,M) struct array with fields x0, relState, C0, tca, cdm, ...
% 
%        .gmmOrder    (1,1) number of Gaussian components per secondary
% 
%        .primary.cart0 (6,1) primary absolute state used for frame conversion
% 
%        .initCovRtn  [bool] true if initial covariances are expressed in RTN
% 
%        .singleObject [bool] duplicate-first-object behavior for multi-object cases
% 
%        .N, .dt, .et  scalars used to generate initial_state.dat for nliGmm
%
% OUTPUT:
%
%   pp = [struct] same input struct with pp.secondary expanded/updated:
% 
%        - each original secondary replaced by gmmOrder entries
% 
%        - for each component: .w, .x0 (ECI) or .relState (RTN), and .C0 updated
%
% EXAMPLE USAGE:
%
%   pp = splitInitial(pp);
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

s = pp.secondary;
n = pp.gmmOrder;

%% Transform into ECI
for k = 1:length(s)
    if pp.singleObject && k > 1
        tca = s(k).tca;
        pp.secondary((k-1)*n+1:k*n) = pp.secondary((k-2)*n+1:(k-1)*n);
        for j = 1:pp.gmmOrder
            pp.secondary((k-1)*n+j).tca = tca;
            pp.secondary((k-1)*n+j).cdm = s(k).cdm;
        end
        continue
    end
    original = s(k);
    if isempty(s(k).x0) && ~isempty(s(k).relState)
        relState = s(k).relState;                                          % [-] (6,1) Relative secondary state in RTN at the time of conjunction
        stateNli = pp.primary.cart0; %%%%%%%%%%%%%%%%%%%%%%%%% da sistemare perché va fatto con il secondario
        [r2ep,wp] = rtn2eci(pp.primary.cart0(1:3,1),pp.primary.cart0(4:6,1));             % [-] (3,3) DCM RTN of primary to ECI for the time of conjunction
        R2Ep      = rot6(r2ep,wp);
        state     = stateNli + R2Ep*relState;
    elseif ~isempty(s(k).x0) && isempty(s(k).relState)
        state = s(k).x0;                                                   % [-] (6,1) Absolute secondary state in ECI at the time of conjunction        
        stateNli  = state;
    else
        error('either the relative or the absolute state of the secondary must be defined');
    end
    [r2es,ws] = rtn2eci(state(1:3),state(4:6));                            % [-] (3,3) DCM RTN of secondary to ECI for the time of conjunction
    R2Es      = rot6(r2es,ws);
    E2Rs      = inv(R2Es);
    if pp.initCovRtn
        C0rtn     = s(k).C0;
        C0        = R2Es*C0rtn*R2Es';
    else
        C0        = s(k).C0;
        C0rtn     = E2Rs*s(k).C0*E2Rs';
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
    [xi, Ci, wi] = VittaldevSplit(xc,C0c,dir,n);
    for j = 1:n
        secInd                  = j + (k-1)*n;           
        pp.secondary(secInd)    = original;
        pp.secondary(secInd).w  = wi(j);
        pp.secondary(secInd).x0 = xi(:,j);                     % state in ECI
        if isempty(s(k).x0) && ~isempty(s(k).relState)
            pp.secondary(secInd).relState = R2Ep\(pp.secondary(secInd).x0 - pp.primary.cart0);
            pp.secondary(secInd).x0 = [];                     % state in ECI
        end
        if pp.initCovRtn
            pp.secondary(secInd).C0 = E2Rs*Ci(:,:,j)*E2Rs';        % covariance in RTN of secondary
        else
            pp.secondary(secInd).C0 = Ci(:,:,j);
        end
    end
end
end