function [A,b_up,b_lo] = altConstr(m,nOpt,alt_const,altMaps,x_exp,pp)
% AltConstr computes linear matrix and the constant term for the altitude
% constraint in MOSEK.
%
% INPUT: m         = [-] (1,1)   Number of optimization variables per node 
%        nOpt      = [-] (1,1)   Total number of optimization variables
%        alt_const = [-] (2,N)   Constant part of the geodetic coordinates
%        altMaps   = [-] (1,6,N) Linear maps linking state to geodetic
%        x_exp     = [-] (6,N)   Expansion point for the linearization
%        pp        = [struct]    Postprocess structure.
% 
% OUTPUT: A    = [-] (2NSKf-2NSK0+1,nOpt) Linear matrix of coefficients of altitude constraint
%         b_up = [-] (2NSKf-2NSK0+1,1)    Upper limit of known terms vector of altitude constraint
%         b_lo = [-] (2NSKf-2NSK0+1,1)    Lower limit of known terms vector of altitude constraint
% 
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------   
%% Define the paramters
ind    = pp.index;                                                         % [struct]  Indeces of the variables per node
n0     = pp.NSK0;                                                          % [-] (1,1) Starting node for SK constraint
nf     = pp.NSKf;                                                          % [-] (1,1) Ending node for SK constraint
nodes  = n0:nf;                                                            % [-] (1,n) Nodes in which the constraint is active
n      = length(nodes);                                                    % [-] (1,1) Number of nodes in which the contraint is active
rows   = 1:n;                                                              % [-] (1,n) Rows in the optimization
altLim = pp.altLim;                                                        % [-] (2,1) Upper and Lower limits for the altitude
altMaps = squeeze(altMaps)';

%% State coefficients
indX1 = reshape(repmat(rows,6,1),6*n,1);                                   % [-] (6n,1) Row indeces to define the state coefficients of the constraints
aux   = repmat((nodes-1)*m,6,1) + ind.state';                              % [-] (6,n)  Auxiliary matrix to construct the column indeces vector
indX2 = reshape(aux,6*n,1);                                                % [-] (6n,1) Column indeces to define the state coefficients of the constraint
valX  = reshape(altMaps(nodes,:)',6*n,1);                                  % [-] (6n,1) Values of the coefficients of the state

%% Virtual buffer coefficients
if pp.enableTrAndVc
    indVc1 = rows';                                                        % [-] (n,1) Row indeces to define the virtual buffer coefficients of the constraints
    aux    = (nodes - 1)*m + ind.altVc';                                   % [-] (1,n)  Auxiliary matrix to construct the column indeces vector
    indVc2 = reshape(aux,n,1);                                             % [-] (n,1) Column indeces to define the virtual buffer coefficients of the constraint
    valVc  = ones(n,1);                                                    % [-] (n,1) Values of the coefficients of the virtual buffer         
else; indVc1 = [];  indVc2 = [];  valVc = [];                              % [-] (0,0)  No virtual control                      
end

%% Build the sparse matrix
ind1 = [indX1; indVc1];                                                    % [-] (7n,1)    Row indeces of full matrix
ind2 = [indX2; indVc2];                                                    % [-] (7n,1)    Column indeces of full matrix
val  = [valX; valVc];                                                      % [-] (7n,1)    Coefficients of full matrix
A    = sparse(ind1,ind2,val);                                              % [-] (6n,nOpt) Sparse matrix of the coefficients
A    = [A sparse(length(rows),nOpt-size(A,2))];                            % [-] (6n,nOpt) Add zeros to have the correct number of columns

%% Build the bounds
expAlt = nan(1,n); 
for i = 1:n; expAlt(:,i) = altMaps(i+n0-1,:)*x_exp(1:6,i+n0-1); end        % [-] (1,n) Transformation of the linearization point 
res   = expAlt - (alt_const(nodes) + 1);                                   % [-] (1,n) Residual of the linerization (subtract one because it is the Earth's radius)
b_up  = res' + altLim(2) + 1;                                              % [-] (n,1) Upper limits
b_lo  = res' + altLim(1) + 1;                                              % [-] (n,1) Lower limits
end