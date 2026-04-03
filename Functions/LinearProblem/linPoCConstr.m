function [A,b_up,b_lo,res] = linPoCConstr(M,m,nOpt,pocMaps,constPoC,x_exp,pp)
% linPoCConstr sets the contraint to define the CA constraint in terms of
% PoC, The  node-wise vector equation that is implemented is:
% 
%              PoC = map*([r_1;...;r_M] -[r_1;...;r_M]_exp) + PoC_const
%  
% INPUT: M           = [-] (1,1) Number of nodes in the problem. 
%        m           = [-] (1,1) Number of optimization variables per node
%        nOpt        = [-] (1,1) Total number of optimization variables.
%        ipcMaps     = [-] (3,M) Maps from relative position to IPC
%        constPc    = [-] (M,1) Constant part of the IPC
%        r_exp       = [-] (3,M) Previous major iteration relative position
%        pp          = [struct]  Postprocess structure
% 
% OUTPUT: A          = [-] (M,nOpt) Linear matrix of coefficients
%         b_up       = [-] (M,1)    Upper limit for the constraint
%         b_lo       = [-] (M,1)    Lower limit for the constraint
%         res        = [-] (1,1)    Residual of the linearization 
% 
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
ind     = pp.index;
indX1   = ones(3*M,1);                                                     % [-] (3M,1) Row indeces to define the relative position coefficients of the constraints
indX2   = nan(3*M,1);
indVc1  = ones(M,1);                                                       % [-] (M,1) Row indeces to define the virtual buffer coefficients of the constraints
indVc2  = nan(M,1);
for j = 1:M
    indX2((j-1)*3 + (1:3)) = (pp.NCA(j) - 1)*m + ind.state(1:3);           % [-] (3M,1) Column indeces to define the relative position coefficients of the constraint
    indVc2(j) = (pp.NCA(j) - 1)*m + ind.caVc(j);                           % [-] (M,1)   Column indeces to define the virtual buffer coefficients of the constraint
end
valX  = reshape(pocMaps',3*M,1);                                           % [-] (3M,1) Values of the coefficients of the relative position
valVc = ones(M,1);                                                         % [-] (M,1)   Values of the coefficients of the virtual buffer

A = sparse([indX1; indVc1], [indX2; indVc2], [valX; valVc]);               % [-] (M,*) Build the sparse matrix                    
A = [A sparse(1,nOpt-size(A,2))];                                          % [-] (M,nOpt) Add zeros to have the correct number of columns

poc_exp   = 0;
for j = 1:M
    poc_exp   = poc_exp + pocMaps(j,:)*x_exp(1:3,pp.NCA(j));                 % [-] (1,1) Expanded part of the PoC
end
res        = poc_exp - constPoC;                                           % [-] (1,1) Residual, comprehensive of the rDeb term to subtract the position of the secondary
b_up       = pp.lim + res;                                                 % [-] (1,1) Upper limit of the constraint
b_lo       = -inf;                                                         % [-] (1,1) Lower limit of the constraint
end