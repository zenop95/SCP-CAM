function [xi, Ci] = propBack(pp)
% PropBack Propagates back the trajectory of the primary s/c to obtain a 
% new state and a new covariance that act as starting point for the 
% optimization. They become the new fixed points for the problem.
% 
% INPUT:  pp = postprocess structure
%       
% OUTPUT: xi = [-] Initial state
%         Ci = [-] initial covariance
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
x0 = pp.primary.cart0;                                                      % [-] (6,1) Starting point to back propagate
N  = pp.N_back+1;                                                          % [-] (1,1) Number of back propagation nodes
et = pp.et + pp.N_back*pp.dt;                                              % [-] (1,1) Ephemeris time at conjucntion

%% Write the .txt file to run the C++ executable
fid = fopen('initial_state.dat', 'w');
fprintf(fid, '%2i\n', N);
fprintf(fid, '%40.12f\n', pp.dt);
fprintf(fid, '%40.12f\n', et);                        
for i = 1:6 
    fprintf(fid, '%40.12f\n', x0(i));
end
fclose(fid);

%% Perform the backward propagation
!wsl ./CppExec/propAidaBack
stm      = dat2maps('maps.dat',6,6,N);                                     % [-] (6,6,N) STMs of the trajectory
x_back   = reshape(load("constPart.dat"),6,N);                             % [-] (6,N)   Back-propagated trajectory
C        = propCovariance(stm,repmat(eye(6),[1,1,N]),N,false,pp);          % [-] (6,6,N) Propagated covariance matrices

%% Initial state and covariance
xi       = pp.cart2x(x_back(:,end));                                       % [-] (6,1) New initial state
Ci       = C(:,:,end);                                                     % [-] (6,6) New initial covariance matrix
end