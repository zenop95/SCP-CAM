function [A,b] = targConstr(m,N,nOpt,target,pp)
% targConstr computes linear matrix and the constant term for the targeting
% constraint in MOSEK.
%
% INPUT: m        = [-] (1,1) Number of optimization variables per node 
%        N        = [-] (1,1) Number of nodes in the optimization
%        nOpt     = [-] (1,1) Total number of optimization
%        target   = [-] (6,1) Relative position final target for the optimization
%        pp       = [struct]  Postprocess structure.
% 
% OUTPUT: A       = [-] (6,nOpt) Linear matrix of coefficients of SK constraint
%         b       = [-] (6,1)    Upper limit of known terms vector of SK constraint
% 
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------   
ind       = pp.index;                                                      % [struct]     Structure of the indeces of each node
indX1     = (1:6)';                                                        % [-] (6,1)    Row indeces to define the state coefficients of the constraint 
indX2     = ind.state' + m*(N-1);                                          % [-] (6,1)    Column indeces to define the state coefficients of the constraint 
valX      = ones(6,1);                                                     % [-] (6,1)    Values of the coefficients of the state
indSlack1 = [(1:6)'; (1:6)'];                                              % [-] (12,1)   Row indeces to define the slack coefficients of the constraint 
indSlack2 = ind.targetSlack';                                              % [-] (12,1)   Column indeces to define the slack coefficients of the constraint 
valSlack  = [ones(6,1); -ones(6,1)];                                       % [-] (6,1)    Values of the coefficients of the slack variables
A         = sparse([indX1; indSlack1],[indX2; indSlack2],[valX; valSlack]);% [-] (6,nOpt) Build the sparse constraint matrix
A         = [A sparse(6,nOpt-size(A,2))];                                  % [-] (6,nOpt) Add zeros to have the correct number of columns
b         = target;                                                        % [-] (6,1)    Known term of the constraint

end