function  ref = refTraj(DAmaps,N,nVars,nDA)
% majorIteration runs the major iteration process for the update of the
% dynamics of the solution.
%
% INPUT: DAmaps = structure of the maps extractd from the DA file
%        N      = number of nodes
%        nVars  = number of variables per node
%        nDA    = number of DA variables
%
% OUTPUT: newTraj   = [m] [m/s] [m/s^2] output relative trajectory
%         majorIter = structure containing output of last minor iteration
%         pp        = structure containing all the simulation input and
%                     output
%         retFlag   = [] (only for MonteCarlo) flag that stops the present 
%                        MonteCarlo sample if the IPC is already below limit 
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

% extract the unperturbed trajectory 
dx0NoMan  = zeros(1,nDA);
ref = nan(nVars,N);

for i = 1:N
    for j = 1:nVars
        % exact position and velocity of the primary at each time step
        ref(j,i) = eval_poly(DAmaps(j,i).C, DAmaps(j,i).E,dx0NoMan);
    end
end