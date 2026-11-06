function pp = GMM_splitInitial(pp)
% GMM_splitInitial Splits each primary & secondary into a Gaussian mixture in ECI/RTN frames
%
% GMM_splitInitial takes the postprocess structure containing the objects
% 
% and replaces each primary & secondary by gmmOrder components following a nonlinearity-
% 
% guided split. For primaries or secondaries given in RTN (relState) the function converts
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
%        .primary     (1,M) struct array with fields x0 (=cart0), C0, tca, cdm, ...
%
%        .secondary   (1,M) struct array with fields x0, relState, C0, tca, cdm, ...
% 
%        .gmmOrder    (1,2) number of Gaussian components per primary and secondary
% 
%        .primary.cart0 (6,1) primary absolute state used for frame conversion
% 
%        .initCovRtn_P and .initCovRtn_S [bool] true if initial covariances are expressed in RTN
% 
%        .singleObject_P and .singleObject_S [bool] duplicate-first-object behavior for multi-object cases
% 
%        .N, .dt, .et  scalars used to generate initial_state.dat for nliGmm
%
% OUTPUT:
%
%   pp = [struct] same input struct with pp.primary and pp.secondary expanded/updated:
% 
%        - each original primary replaced by gmmOrder(1) entries
%
%        - each original secondary replaced by gmmOrder(2) entries
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

p = pp.primary;
n1 = pp.gmmOrder(1);

s = pp.secondary;
n2 = pp.gmmOrder(2);

%% Transform primary into ECI
for k = 1:length(p)
    %if pp.singleObject_P && k > 1
    if k > 1
        tca = p(k).tca;
        pp.primary((k-1)*n1+1:k*n1) = pp.primary((k-2)*n1+1:(k-1)*n1);
        for j = 1:pp.gmmOrder(1)
            pp.primary((k-1)*n1+j).tca = tca;
            pp.primary((k-1)*n1+j).cdm = p(k).cdm;
        end
        continue
    end
    original = p(k);
    if ~isempty(p(k).x0)   % && isempty(p(k).relState)
        state = p(k).x0;                                                  % [-] (6,1) Absolute primary state in ECI at the time of conjunction        
        stateNli  = state;
    else
        error('absolute state of the primary must be defined');
    end
    % if isempty(s(k).x0) && ~isempty(s(k).relState)
    %     relState = s(k).relState;                                          % [-] (6,1) Relative secondary state in RTN at the time of conjunction
    %     stateNli = pp.primary.cart0; %%%%%%%%%%%%%%%%%%%%%%%%% da sistemare perché va fatto con il secondario
    %     [r2ep,wp] = rtn2eci(pp.primary.cart0(1:3,1),pp.primary.cart0(4:6,1));             % [-] (3,3) DCM RTN of primary to ECI for the time of conjunction
    %     R2Ep      = rot6(r2ep,wp);
    %     state     = stateNli + R2Ep*relState;
    % elseif ~isempty(s(k).x0) && isempty(s(k).relState)
    %     state = s(k).x0;                                                   % [-] (6,1) Absolute secondary state in ECI at the time of conjunction        
    %     stateNli  = state;
    % else
    %     error('either the relative or the absolute state of the secondary must be defined');
    % end
    [r2ep,wp] = rtn2eci(state(1:3),state(4:6));                            % [-] (3,3) DCM RTN of secondary to ECI for the time of conjunction
    R2Ep      = rot6(r2ep,wp);
    E2Rp      = inv(R2Ep);
    if pp.initCovRtn_P
        C0rtn     = p(k).C0;
        C0        = R2Ep*C0rtn*R2Ep';
    else
        C0        = p(k).C0;
        C0rtn     = E2Rp*p(k).C0*E2Rp';
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
    [xi, Ci, wi] = VittaldevSplit_Multipledir(xc,C0c,dir,n1);
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
end


%% Transform secondary into ECI
for k = 1:length(s)
    %if pp.singleObject_S && k > 1
    if k > 1
        tca = s(k).tca;
        pp.secondary((k-1)*n2+1:k*n2) = pp.secondary((k-2)*n2+1:(k-1)*n2);
        for j = 1:pp.gmmOrder(2)
            pp.secondary((k-1)*n2+j).tca = tca;
            pp.secondary((k-1)*n2+j).cdm = s(k).cdm;
        end
        continue
    end
    original = s(k);
    if isempty(s(k).x0) && ~isempty(s(k).relState)
        relState = s(k).relState;                                          % [-] (6,1) Relative secondary state in RTN at the time of conjunction
        stateNli = pp.primaryReferenceState; %%%%%%%%%%%%%%%%%%%%%%%%% da sistemare perché va fatto con il secondario
        [r2ep,wp] = rtn2eci(stateNli(1:3,1),stateNli(4:6,1));             % [-] (3,3) DCM RTN of primary to ECI for the time of conjunction
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
    if pp.initCovRtn_S
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
end