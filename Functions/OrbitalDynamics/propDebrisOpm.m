function pp = propDebrisOpm(pp,x0,N_forw,N_back)
% propDebris Propagates the debris orbit for the required time according to
% one of the possible dynamics models. For the J2 and Full dynamics model
% it needs to call the corresponding c++ executable.
% INPUT:
%       pp: parameters straucture
%       x0: [m][m/s] Initial state of the debris.
%       t:  [s] time history
%
% OUTPUT:
%       x_dProp: [m][m/s] evolution of the state of the debris
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
% write initial conditions and variables to file to pass to C++
%% Forward propagation
N = N_forw+1;
fid = fopen('initial_state.dat', 'w');
fprintf(fid, '%2i\n', N);
fprintf(fid, '%40.12f\n', pp.dt);
fprintf(fid, '%40.12f\n', pp.et + N_back*pp.dt);
fprintf(fid, '%2i\n', 1); %forward flag
for i = 1:6 
    fprintf(fid, '%40.12f\n', x0(i));
end
fclose(fid);
!wsl ./CppExec/propAidaDebris
pp.secondary.x = reshape(load('constPart.dat'),6,N);
% Polynomials organized as state vectors
mapsDebForw    = dat2maps('maps.dat',6,6,N);
covForw        = nan(6,6,N);
covForw(:,:,1) = pp.C0s;
for i = 2:N
    covForw(:,:,i) = mapsDebForw(:,:,i)*covForw(:,:,i-1)*mapsDebForw(:,:,i)';
end
mapsDebBack = nan(6,6,0);
covBack     = nan(6,6,0);
%% Backward propagation
if N_back > 0
    N = N_back+1;
    fid = fopen('initial_state.dat', 'w');
    fprintf(fid, '%2i\n',     N);
    fprintf(fid, '%40.12f\n', pp.dt);
    fprintf(fid, '%40.12f\n', pp.et + N_back*pp.dt);
    fprintf(fid, '%2i\n', 2); %backward flag
    for i = 1:6 
        fprintf(fid, '%40.12f\n', x0(i));
    end
    fclose(fid);
    !wsl ./CppExec/propAidaDebris
    ref = reshape(load("constPart.dat"),6,N);
    % Polynomials organized as state vectors
    x_dPropBack     = flip(ref(:,2:end),2); % exclude first node redundant with forw prop
    pp.secondary.x  = [x_dPropBack, pp.secondary.x];
    pp.secondary.x0 = pp.secondary.x(:,1);
    mapsDebBack     = dat2maps('maps.dat',6,6,N);
    covBack         = nan(6,6,N);
    covBack(:,:,1)  = pp.C0s;
    for i = 2:N
        covBack(:,:,i) = mapsDebBack(:,:,i)*covBack(:,:,i-1)*mapsDebBack(:,:,i)'; 
    end
end
pp.secondary.stm        = cat(3,flip(permute(mapsDebBack(:,:,2:end),...
                                                  [2 1 3]),3),mapsDebForw); % transpose and flip the stm to have it in the forward direction
pp.secondary.covariance = cat(3,flip(covBack(:,:,2:end),3),covForw);        % flip the covariance to have it in the forward direction
end