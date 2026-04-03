function [A,b_up,b_lo] = nliConstr(N,m,nOpt,x_exp,nu,ind)
% NliConstr sets the contraint to define the NLI variables used to
% impose the trust region. The  node-wise and component-wise equation that  
% is implemented is:
% 
%                   nli = (x - x_exp)*nu
%
% INPUT: N           = [-] (1,1) Number of nodes in the problem. 
%        m           = [-] (1,1) Number of optimization variables per node
%        nOpt        = [-] (1,1) Total number of optimization variables.
%        x_exp       = [-] (9,N) Expansion point for the linearization
%        nu          = [-] (6,N) Non linear parameter
%        ind         = [struct] Structure that contains the indexes
% 
% OUTPUT: A          = [-] (6N,nOpt) Linear matrix of coefficients
%         b_up       = [-] (6N,1)    Upper limits for the constraints
%         b_lo       = [-] (6N,1)    Lower limits for the constraints
% 
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
rows    = 1:6*N;                                                           % [-] (1,6N) Number of rows for the dynamics A matrix
indX1   = rows';                                                           % [-] (6N,1) Row indices to define the state coefficients of the constraints
aux     = repmat((0:1:N-1)*m,6,1) + ind.state';                            % [-] (6,N)  Auxiliary array to define the column indices
indX2   = reshape(aux,6*N,1);                                              % [-] (6N,1) Column indices to define the state coefficients of the constraints 
valX    = reshape(nu(1:6,:),6*N,1);                                        % [-] (6N,1) Values of the state coefficients of the contraints
indNli1 = rows';                                                           % [-] (6N,1) Row indices to define the NLI coefficients of the constraints
aux     = repmat((0:1:N-1)*m,6,1) + ind.tr';                               % [-] (6,N)  Auxiliary array to define the column indices 
indNli2 = reshape(aux,6*N,1);                                              % [-] (6N,1) Column indices to define the NLI coefficients of the constraints
valNli  = ones(6*N,1);                                                     % [-] (6N,1) Values of the NLI coefficients of the contraints 
A       = sparse([indX1; indNli1],[indX2; indNli2],[valX; valNli]);        % [-] (6N,*) Build the sparse matrix of the coefficients
A       = [A sparse(length(rows),nOpt-size(A,2))];                         % [-] (6N,nOpt) Add zeros to have the correct number of columns
b_up    = reshape(x_exp(1:6,:).*nu(1:6,:),6*N,1);                          % [-] (6N,1) Create the upper limits vector of the contraints
b_lo    = b_up;                                                            % [-] (6N,1) Create the equality constraint setting the lower equal to the upper limit
end