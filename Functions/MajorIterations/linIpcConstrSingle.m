function [A,b_up,b_lo] = linIpcConstrSingle(N,M,m,nOpt,ipcMaps,constIpc,x_exp,pp)
% linIpcConstr sets the contraint to define the CA constraint in terms of
% IPC, The  node-wise vector equation that is implemented is:
% 
%                   IPC < map*(p - p_exp) + IPC_const
%  
% INPUT: N           = [-] (1,1) Number of nodes in the problem. 
%        m           = [-] (1,1) Number of optimization variables per node
%        nOpt        = [-] (1,1) Total number of optimization variables.
%        ipcMaps     = [-] (3,N) Maps from relative position to IPC
%        constIpc    = [-] (N,1) Constant part of the IPC
%        r_exp       = [-] (3,N) Previous major iteration relative position
%        pp          = [struct]  Postprocess structure
% 
% OUTPUT: A          = [-] (N,nOpt) Linear matrix of coefficients
%         b_up       = [-] (N,1)    Upper limit for the constraint
%         b_lo       = [-] (N,1)    Lower limit for the constraint
% 
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
[a,b]   = max(constIpc);
rows    = 1;                                                               % [-] (1,N)   Number of scalar constraints
ind     = pp.index;
indPos1 = ones(3,1);                                                       % [-] (3,1) Row indeces to define the relative position coefficients of the constraints
indVc1  = 1;                                                               % [-] (1,1) Row indeces to define the virtual buffer coefficients of the constraints
indPos2 = (b - 1)*m + ind.state(1:3)';                                     % [-] (3,1) Column indeces to define the relative position coefficients of the constraint
indVc2  = (b - 1)*m + ind.caVc(1);                                         % [-] (1,1)   Column indeces to define the virtual buffer coefficients of the constraint
valPos  = sum(ipcMaps(:,:,b),1)';                                          % [-] (3,1) Values of the coefficients of the relative position
valVc   = 1;                                                               % [-] (1,1)   Values of the coefficients of the virtual buffer

A = sparse([indPos1; indVc1], [indPos2; indVc2], [valPos; valVc]);% [-] (N,*)    Build the sparse matrix                    
A = [A sparse(length(rows),nOpt-size(A,2))];                               % [-] (N,nOpt) Add zeros to have the correct number of columns

ipc_exp  = sum(ipcMaps(:,:,b),1)*x_exp(1:3,b);                             % [-] (1,1) Expanded part of the TIPC
b_up     = pp.lim + ipc_exp - a;                                           % [-] (N,1) Residual of the linearization for upper limit
b_lo     = -inf;                                                           % [-] (N,1) Lower limit of the constraints
end