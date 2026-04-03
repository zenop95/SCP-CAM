function [secondary,NCA] = propSecondary(secondary,state,NCA,pp)
% PropSecondary Propagates the secondary orbit for the required time 
% and it propagates the s/c covariance with the associatedd linear maps. 
% 
% INPUT:  x0 = [-] (6,1) State of the secondary spacecraft at conjunction
%         pp = [struct]  Postprocess structure
%
% OUTPUT: pp = [struct]  Postprocess structure
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
Ntca       = secondary.tca;
cdm        = secondary.cdm;
NCA        = [NCA; Ntca + pp.N_back];                                      % [-] (1,1) Conjunction nodes
n          = length(NCA); 
if pp.singleObject && n > pp.gmmOrder
    secondary     = pp.secondary(n-pp.gmmOrder);
    secondary.tca = Ntca;
    secondary.cdm = cdm;
    return
end
[r2e,w]    = rtn2eci(state(1:3,NCA(end)),state(4:6,NCA(end)));             % [-] (3,3) DCM RTN of primary to ECI for the time of conjunction
R2E        = rot6(r2e,w);
if isempty(secondary.x0) && ~isempty(secondary.relState) && ~secondary.ang
    x0 = R2E*secondary.relState + state(:,NCA(end));                       % [-] (6,1) Secondary state in ECI at the time of conjunction
    secondary.x0 = x0;
elseif isempty(secondary.x0) && ~isempty(secondary.relState) && secondary.ang
    rotAx = r2e*[1;0;0];
    x0 = [r2e*secondary.relState(1:3) + state(1:3,NCA(end)); 
          axang2dcm(rotAx,secondary.velAng)*state(4:6,NCA(end))];          % [-] (6,1) Secondary state in ECI at the time of conjunction
    secondary.x0 = x0;
elseif ~isempty(secondary.x0)
    x0 = secondary.x0;                                                     % [-] (6,1) Secondary state in ECI at the time of conjunction
else
    error('either the relative or the absolute state of the secondary must be defined');
end

if pp.initCovRtn 
    [r2es,ws] = rtn2eci(x0(1:3),x0(4:6));                                      % [-] (3,3) DCM RTN of secondary to ECI for the time of conjunction
    R2Es      = rot6(r2es,ws);
    C0        = R2Es*secondary.C0*R2Es';
else
    C0        = secondary.C0;
end
%% Forward propagation
[x, cov] = mapsSecondary(Ntca,x0,C0,pp);                                   % [-] (6,N) (6,6,N) Propagate forward the trajectory and the covariance

%% Covariancce as output
secondary.x          = x;                                                  % [-] (6,N)   Full trajectory
secondary.covariance = cov;                                                % [-] (6,6,N) Full covariance matrices
end