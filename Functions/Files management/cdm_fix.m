function cdm = cdm_fix(pp,cdm,n)
% Cdm_fix writes a structure where the relevant parameters extracted from
% the Conjunction Data Message (CDM) are stored
%
% INPUT: cdm = CDM structure from Andrea's function (read_cdm)
% 
%        n   = Node in which the conjunction happens
%
% OUTPUT: cdm = New CDM structure
%
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

Lsc = 6378.173;                                                            % [m] distance scaling factor
Vsc = sqrt(398600.4418/Lsc);                                               % [m/s] velocity scaling factor
scaling = [Lsc*ones(3,1); Vsc*ones(3,1)];                                  % vector with scaling parameters to scale state vectors
D = diag(1./scaling);                                                      % Diagonal matrix to scale covariances
tca_spl      = split(cdm.TCA,{'T','-',':'})';
cdm.tca_spl  = tca_spl([3,2,1,4,5,6]);                                     % [utc] UTC of the conjunction
cdm.et       = pp.opm.et(n);                                               % [s] ephemeris time of the conjunction
cdm.relState = (cdm.object1.data.mean_state - ...
                                cdm.object2.data.mean_state)./scaling;     % [-] Relative state at conjunction
state1       = pp.opm.nodeStates(:,n);                                     % [-] state of the primary s/c at conjunction
state2       = state1-cdm.relState;                                        % [-] state of the secondary s/c at conjunction
cov1Rtn      = D*cdm.object1.data.covariance*D;                            % [-] primary s/c covariance matrix 
cov2Rtn      = D*cdm.object2.data.covariance*D;                            % [-] secondary s/c covariance matrix 

%% Compute transform matrices
[lvlh2eci1,w1] = lvlh2eci(state1(1:3),state1(4:6));                        % DCM and angular velocity between ECI and LVLH for primary s/c
rtn2eci1       = lvlh2eci1*lvlh2rtn';
wRtn2eci1      = lvlh2rtn*w1;
% transf1        = [rtn2eci1 zeros(3,3); skew(wRtn2eci1) rtn2eci1];
transf1        = [rtn2eci1 zeros(3,3); zeros(3,3) rtn2eci1];               % Transform matrix between RTN and ECI for primary s/c

[lvlh2eci2,w2] = lvlh2eci(state2(1:3),state2(4:6));                        % DCM and angular velocity between ECI and LVLH for primary s/c
rtn2eci2       = lvlh2eci2*lvlh2rtn';
wRtn2eci2      = lvlh2rtn*w2;
% transf2        = [rtn2eci2 zeros(3,3); skew(wRtn2eci2) rtn2eci2];
transf2        = [rtn2eci2 zeros(3,3); zeros(3,3) rtn2eci2];               % Transform matrix between RTN and ECI for secondary s/c

%% everything in ECI 
cdm.cov1  = transf1*cov1Rtn*transf1';                                      % [-] ECI covariance of primary s/c at conjunction
cdm.cov2  = transf2*cov2Rtn*transf2';                                      % [-] ECI covariance of secondary s/c at conjunction


end

