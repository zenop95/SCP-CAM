function secondary = propSecondary(secondary, pp)
% PropSecondary Propagates the secondary's orbit for the required time 
% and it propagates the s/c covariance with the associated linear maps. 
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
x0 = secondary.x0(:);
C0 = secondary.C0;

%% Forward propagation
[x, cov] = mapsSecondary(Ntca,x0,C0,pp);                                   % [-] (6,N) (6,6,N) Propagate forward the trajectory and the covariance

%% Covariancce as output
secondary.x          = x;                                                  % [-] (6,N)   Full trajectory
secondary.cart       = x;
secondary.covariance = cov;                                                % [-] (6,6,N) Full covariance matrices
end