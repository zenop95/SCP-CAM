function [A,b_up,b_lo] = skConstr(m,nOpt,skLim,ll_const,llMaps,x_exp,pp)
% SkConstr computes linear matrix and the constant term for the SK
% constraint in MOSEK. This also includes the final target constraint 
% depending on what is present in the pp structure.
%
% INPUT: m        = [-] (1,1)   Number of optimization variables per node 
%        nOpt     = [-] (1,1)   Total number of optimization variables
%        skLim    = [-] (2,1)   Vector with amplitude of the SK box in
%        lla      = [-] (2,N)   Constant part of the geodetic coordinates
%        llMaps   = [-] (2,6,N) Linear maps linking state to geodetic
%        x_exp    = [-] (6,N)   Expansion point for the linearization
%        pp       = [struct]    Postprocess structure.
% 
% OUTPUT: A    = [-] (2NSKf-2NSK0+1,nOpt) Linear matrix of coefficients of SK constraint
%         b_up = [-] (2NSKf-2NSK0+1,1)    Upper limit of known terms vector of SK constraint
%         b_lo = [-] (2NSKf-2NSK0+1,1)    Lower limit of known terms vector of SK constraint
% 
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------   
%% Define the paramters
ind    = pp.index;                                                         % [struct]   Indeces of the variables per node
n0     = pp.NSK0;                                                          % [-] (1,1)  Starting node for SK constraint
nf     = pp.NSKf;                                                          % [-] (1,1)  Ending node for SK constraint
nodes  = n0:nf;                                                            % [-] (1,n)  Nodes in which the constraint is active
n      = length(nodes);                                                    % [-] (1,1)  Number of nodes in which the contraint is active
lla0   = [0; pp.nomLon];                                                   % [-] (2,1)  Nominal value of the geodetic coordinates
rows   = 1:2*n;                                                            % [-] (1,2n) Rows in the optimization

%% State coefficients
indX1 = reshape(repmat(rows,6,1),12*n,1);                                  % [-] (12n,1) Row indeces to define the state coefficients of the constraints
aux   = repmat(repmat((nodes-1)*m,6,1) + ind.state',2,1);                  % [-] (12,n)  Auxiliary matrix to construct the column indeces vector
indX2 = reshape(aux,12*n,1);                                               % [-] (12n,1) Column indeces to define the state coefficients of the constraint
valX  = reshape(permute(llMaps(:,:,nodes),[2 1 3]),12*n,1);                % [-] (12n,1) Values of the coefficients of the state

%% Virtual buffer coefficients
indVc1 = rows';                                                        % [-] (2n,1) Row indeces to define the virtual buffer coefficients of the constraints
aux    = repmat((nodes-1)*m,2,1) + ind.skVc';                          % [-] (2,n)  Auxiliary matrix to construct the column indeces vector
indVc2 = reshape(aux,2*n,1);                                           % [-] (2n,1) Column indeces to define the virtual buffer coefficients of the constraint
valVc  = ones(2*n,1);                                                  % [-] (2n,1) Values of the coefficients of the virtual buffer         


%% Slack variables coefficients
if pp.skSoft
    indSlack1 = reshape(repmat(rows,4,1),8*n,1);                           % [-] (8n,1) Row indeces to define the slack variables coefficients of the constraints
    aux       = repmat(repmat((nodes-1)*m,4,1) + ind.skSlack',2,1);        % [-] (8,n)  Auxiliary matrix to construct the column indeces vector
    indSlack2 = reshape(aux,8*n,1);                                        % [-] (8n,1) Column indeces to define the slack variables  coefficients of the constraint
    valSlack  = ones(8*n,1);                                               % [-] (8n,1) Values of the coefficients of the slack variables        
else; indSlack1 = [];  indSlack2 = [];  valSlack = [];                     % [-] (0,0)  No slack variables
end

%% Build the sparse matrix
ind1 = [indX1; indVc1; indSlack1];                                         % [-] (22n,1)    Row indeces of full matrix
ind2 = [indX2; indVc2; indSlack2];                                         % [-] (22n,1)    Column indeces of full matrix
val  = [valX; valVc; valSlack];                                            % [-] (22n,1)    Coefficients of full matrix
A    = sparse(ind1,ind2,val);                                              % [-] (12n,nOpt) Sparse matrix of the coefficients
A    = [A sparse(length(rows),nOpt-size(A,2))];                            % [-] (12n,nOpt) Add zeros to have the correct number of columns

%% Build the bounds
expLatLon = nan(2,n); 
for i = 1:n; expLatLon(:,i) = llMaps(:,:,i+n0-1)*x_exp(1:6,i+n0-1); end    % [-] (2,n)   Transformation of the linearization point
res   = expLatLon - ll_const(:,nodes);                                     % [-] (2,n)   Residual of the linerization
limUp = (lla0 + skLim).*ones(2,n);                                         % [-] (2,n)   Upper boundary of the SK box
limLo = (lla0 - skLim).*ones(2,n);                                         % [-] (2,n)   Lower boundary of the SK box
b_up  = reshape(res + limUp,2*n,1);                                        % [-] (12n,1) Upper limits
b_lo  = reshape(res + limLo,2*n,1);                                        % [-] (12n,1) Lower limits
end