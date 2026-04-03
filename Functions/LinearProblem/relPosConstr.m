function [A,b] = relPosConstr(N,m,conv,nOpt,constPos,x_exp,j,pp)
% relPosConstr sets the contraint to define the relative position variables
% in cartesian coordinates. The  node-wise vector equation that is 
% implemented is:
% 
%                   dr = B_m2c*(x - x_exp) + r_p - r_s
%
%  
% INPUT: N           = [-] (1,1)   Number of nodes in the problem. 
%        m           = [-] (1,1)   Number of optimization variables per node
%        conv        = [-] (6,6,N) Conversion maps from MEE to cartesian
%        nOpt        = [-] (1,1)   Total number of optimization variables.
%        rDeb        = [-] (3,N)   Cartesian position of the secondary
%        constPos    = [-] (3,N)   Constant part of th transformation from state to cartesian
%        x_exp       = [-] (6,N)   Previous major iteration solution
%        pp          = [struct]    Postprocess structure
% 
% OUTPUT: A          = [-] (3N,nOpt) Linear matrix of coefficients
%         b_up       = [-] (3N,1)    Equality limit for the constraint
% 
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
ind        = pp.index;
rDeb       = pp.secondary(j).cart(1:3,:);
rows       = 1:3*N;                                                        % [-] (1,3N) Number of scalar constraints
indX1      = reshape(repmat(rows,6,1),18*N,1);                             % [-] (3N,1) Row indeces to define the state coefficients of the constraints
aux        = repmat(repmat((0:1:N-1)*m,6,1) + ind.state',3,1);             % [-] (3,N)  Auxiliary matrix to construct the column indeces vector
indX2      = reshape(aux,18*N,1);                                          % [-] (3,N)  Column indeces to define the state coefficients of the constraint
valX       = reshape(permute(conv(1:3,:,:),[2 1 3]),18*N,1);               % [-] (3,N)  Values of the coefficients of the state (maps that convert from state to cartesian)
indRelPos1 = rows';                                                        % [-] (3,N)  Row indeces to define the relative position coefficients of the constraint
aux        = repmat((0:1:N-1)*m,3,1) + ind.relPos((1:3)+(j-1)*3)';         % [-] (3N,1) Auxiliary matrix to construct the column indeces vector
indRelPos2 = reshape(aux,3*N,1);                                           % [-] (3,N)  Column indeces to define the relative position coefficients of the constraint
valRelPos  = -ones(3*N,1);                                                 % [-] (3,N)  Values of the coefficients of the relative position
A          = sparse([indX1; indRelPos1],[indX2; indRelPos2], ...           % [-] (6N,*) Build the sparse matrix
                    [valX; valRelPos]);
A          = [A sparse(length(rows),nOpt-size(A,2))];                      % [-] (6N,nOpt) Add zeros to have the correct number of columns
x_exp      = [reshape([x_exp; zeros(m-6,N)],m*N,1); zeros(pp.sl,1)];       % [-] (6N,*)    Expansion point vector
cart_exp   = A*x_exp;                                                      % [-] (6N,1)    Expanded part of the residual of the linearization
b          = cart_exp + reshape(rDeb - constPos,3*N,1);                    % [-] (6N,1)    Residual, comprehensive of the rDeb term to subtract the position of the secondary
end