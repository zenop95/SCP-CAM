function [A, b] = contConstrJit(N,m,sl,dynMaps,x_const,x_exp,pp)
% smdGrdContr computes linear matrix and the constant term for the SMD
% gradient constraint in MOSEK
%
% INPUT: N       = [-] (1x1)   Number of nodes in the problem. 
%        m       = [-] (1x1)   Number of optimization variables per node
%        dynMaps = [-] (6x9xN) Matrices of the linear maps between nodes
%        x_const = [-] (6xN)   Constant part of the propagation
%        x_exp   = [-] (6xN)   Expansion point for the linearization
%        pp      = [struct]    Postprocess structure
% 
% OUTPUT: A      = [-] (6N,nOpt) Linear matrix of coefficients
%         b_up   = [-] (6N,1)    Known terms of the constraints
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------   
ind      = pp.index;                                                       % [struct]   Indeces of the variables per node
nOpt     = m*N + sl;                                                       % [-] (1,1)  Number of optimization variables
n        = N-1;                                                            % [-] (1,1)  Number of nodes where to apply the constraint
rows     = 1:6*n;                                                          % [-] (1,6n) Rows in the optimization

%% State coefficients
indX1    = reshape(repmat(rows,6,1),36*n,1);                               % [-] (36n,1) Row indeces to define the state coefficients of the constraints
aux      = repmat(repmat((0:1:n-1)*m,6,1) + ind.state',6,1);               % [-] (36,n)  Auxiliary matrix to construct the column indeces vector
indX2    = reshape(aux,36*n,1);                                            % [-] (36n,1) Column indeces to define the state coefficients of the constraint
valX     = reshape(permute(dynMaps(:,1:6,2:end),[2 1 3]),36*n,1);          % [-] (36n,1) Values of the coefficients of the state

%% Throttle variable coefficients
Bu = nan(6,n);
for i = 1:n
    Bu(:,i) = dynMaps(:,7:9,i+1)*pp.uFix;
end
indTh1    = reshape(rows,6*n,1);                                           % [-] (6n,1) Row indeces to define the throttle coefficients of the constraints
aux       = repmat((0:1:n-1)*m,6,1) + ind.ctrl;                            % [-] (6,n)  Auxiliary matrix to construct the column indeces vector
indTh2    = reshape(aux,6*n,1);                                            % [-] (6n,1) Column indeces to define the throttle coefficients of the constraint
valTh     = reshape(Bu,6*n,1);                                             % [-] (6n,1) Values of the coefficients of the throttle


%% Propagated state coefficients
indNext1 = rows';                                                          % [-] (6n,1) Row indeces to define the propagate state coefficients of the constraints
aux      = repmat((1:1:n)*m,6,1) + ind.state';                             % [-] (6,n)  Auxiliary matrix to construct the column indeces vector
indNext2 = reshape(aux,6*n,1);                                             % [-] (6n,1) Column indeces to define the propagate state coefficients of the constraint
valNext  = -ones(length(indNext1),1);                                      % [-] (6n,1) Values of the coefficients of the propagate state

%% Virtual control coefficients
indVc1 = rows';                                                        % [-] (6n,1) Row indeces to define the virtual control coefficients of the constraints
aux    = repmat((1:1:n)*m,6,1) + ind.stateVc';                         % [-] (6,n)  Auxiliary matrix to construct the column indeces vector
indVc2 = reshape(aux,6*n,1);                                           % [-] (6n,1) Column indeces to define the virtual control coefficients of the constraint
valVc  = ones(length(indVc1),1);                                       % [-] (6n,1) Values of the coefficients of the virtual control


%% Build constraint matrix and limits
ind1       = [indX1; indTh1; indNext1; indVc1];                            % [-] (66,n) Full row indeces for the constraint matrix
ind2       = [indX2; indTh2; indNext2; indVc2];                            % [-] (66,n) Full column indeces for the constraint matrix
val        = [valX;  valTh;  valNext;  valVc];                             % [-] (66,n) Values of the coefficients of the constraint matrix
A          = sparse(ind1,ind2,val);                                        % [-] (6n,*) Build the sparse constraint matrix
A          = [A sparse(length(rows),nOpt-size(A,2))];                      % [-] (6n,nOpt) Add zeros to have the correct number of columns
A1         = sparse([indX1; indTh1],[indX2; indTh2],[valX; valTh]);        % [-] (6n,*) Build the marix for the propagation of the expansion point
A1         = [A1 sparse(length(rows),nOpt-size(A1,2))];                    % [-] (6n,nOpt) Add zeros to have the correct number of columns
x_exp      = [reshape([x_exp; zeros(m-7,N)],m*N,1); zeros(sl,1)];          % [-] (mN,1) Reshape expansion point for propagation
prop_exp   = A1*x_exp;                                                     % [-] (6n,1) Propagation of the expanded point
const_part = reshape(x_const(:,2:end),6*n,1);                              % [-] (6n,1) Reshape constant part of the expansion
b          = prop_exp - const_part;                                        % [-] (6n,1) Residual of the linearization
end