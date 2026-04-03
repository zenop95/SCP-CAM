function [A,b] = relPosConstrAvA(N,m,nOpt,p,x_exp,pp)
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
conv1     = p(1).cartMaps;
conv2     = p(2).cartMaps;
constPos1 = p(1).cart(1:3,:);
constPos2 = p(2).cart(1:3,:);

ind1       = pp.index(1);
ind2       = pp.index(2);
rows       = 1:3*N;                                                        % [-] (1,3N) Number of scalar constraints
indX11     = reshape(repmat(rows,6,1),18*N,1);                             % [-] (3N,1) Row indeces to define the state coefficients of the constraints
indX12     = indX11;                                                       % [-] (3N,1) Row indeces to define the state coefficients of the constraints
aux        = repmat(repmat((0:1:N-1)*m,6,1) + ind1.state',3,1);            % [-] (3,N)  Auxiliary matrix to construct the column indeces vector
indX21     = reshape(aux,18*N,1);                                          % [-] (3,N)  Column indeces to define the state coefficients of the constraint
aux        = repmat(repmat((0:1:N-1)*m,6,1) + ind2.state',3,1);            % [-] (3,N)  Auxiliary matrix to construct the column indeces vector
indX22     = reshape(aux,18*N,1);                                          % [-] (3,N)  Column indeces to define the state coefficients of the constraint
valX1      = reshape(permute(conv1(1:3,:,:),[2 1 3]),18*N,1);             % [-] (3,N)  Values of the coefficients of the state (maps that convert from state to cartesian)
valX2      = -reshape(permute(conv2(1:3,:,:),[2 1 3]),18*N,1);              % [-] (3,N)  Values of the coefficients of the state (maps that convert from state to cartesian)
indRelPos1 = rows';                                                        % [-] (3,N)  Row indeces to define the relative position coefficients of the constraint
aux        = repmat((0:1:N-1)*m,3,1) + ind1.relPos';                       % [-] (3N,1) Auxiliary matrix to construct the column indeces vector
indRelPos2 = reshape(aux,3*N,1);                                           % [-] (3,N)  Column indeces to define the relative position coefficients of the constraint
valRelPos  = -ones(3*N,1);                                                 % [-] (3,N)  Values of the coefficients of the relative position
A          = sparse([indX11; indX12; indRelPos1], ...
                [indX21; indX22; indRelPos2],[valX1; valX2; valRelPos]);   % [-] (6N,*) Build the sparse matrix
                    
A          = [A sparse(length(rows),nOpt-size(A,2))];                      % [-] (6N,nOpt) Add zeros to have the correct number of columns
x_exp     = [reshape([x_exp(:,:,1); ...
                zeros(m/2-(6+length(ind1.ctrl)),N); ...
                x_exp(:,:,2); zeros(m/2-(6+length(ind2.ctrl)),N)] ...
                ,m*N,1); zeros(pp.sl,1)];                                  % [-] (6N,*)    Expansion point vector
cart_exp  = A*x_exp;                                                       % [-] (6N,1)    Expanded part of the residual of the linearization
b         = cart_exp + reshape(constPos2 - constPos1,3*N,1);               % [-] (6N,1)    Residual, comprehensive of the rDeb term to subtract the position of the secondary
end