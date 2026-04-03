function x_i = backPropagate(pp,x0,tf)
% backPropagate.
% INPUT:
%        x_0 = initial state
%        t   = back propagation time
% OUTPUT:
%        x_i: [m][m/s] back propagated state of the debris
%
% Reference: Malyuta, D., & Acikmese, B. (2021). Fast Homotopy for 
% Spacecraft Rendezvous Trajectory Optimization with Discrete Logic. 
% http://arxiv.org/abs/2107.07001
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

% write initial conditions and variables to file to pass to C++
fid = fopen('initial_state.dat', 'w');
    for i = 1:6 
        fprintf(fid, '%40.12f\n', x0(i));
    end
    fprintf(fid, '%40.12f\n', tf);
    fprintf(fid, '%40.12f\n', pp.dt);
    fprintf(fid, '%40.12f\n', pp.jd);
fclose(fid);

!wsl ./CppExec/propAidaBack
x_i  = load('maps.dat');
end