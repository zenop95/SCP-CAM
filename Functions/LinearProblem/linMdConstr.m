function [A,b_up,b_lo] = linMdConstr(M,m,nOpt,mdMaps,posConst,r_exp,pp)
% linMdConstr sets the contraint to define the CA constraint in terms of
% miss distance. 
%  
% INPUT: M           = [-] (1,1) Number of nodes in the problem. 
%        m           = [-] (1,1) Number of optimization variables per node
%        nOpt        = [-] (1,1) Total number of optimization variables.
%        mdMaps      = [-] (3,M) Maps from relative position to miss
%                                distance
%        posConst    = [-] (M,1) Constant part of the absolute position
%        r_exp       = [-] (3,M) Previous major iteration relative position
%        pp          = [struct]  Postprocess structure
% 
% OUTPUT: A          = [-] (M,nOpt) Linear matrix of coefficients
%         b_up       = [-] (M,1)    Upper limit for the constraint
%         b_lo       = [-] (M,1)    Lower limit for the constraint
% 
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
ind        = pp.index;
indRelPos1 = repmat(1:M,[3,1]);                                            % [-] (3M,1) Row indeces to define the relative position coefficients of the constraints
indRelPos2 = nan(3*M,1);
indVc1     = ones(M,1);                                                    % [-] (M,1) Row indeces to define the virtual buffer coefficients of the constraints
indVc2     = nan(M,1);
md_exp     = nan(M,1);
mdConst    = nan(M,1);
for j = 1:M
    indRelPos2((j-1)*3 + (1:3)) = (pp.NCA(j) - 1)*m + ...
                                        ind.relPos((j-1)*3 + (1:3));       % [-] (3M,1) Column indeces to define the relative position coefficients of the constraint
    indVc2(j) = (pp.NCA(j) - 1)*m + ind.caVc(j);                           % [-] (M,1)   Column indeces to define the virtual buffer coefficients of the constraint
    md_exp(j) = mdMaps(j,:)*r_exp(:,pp.NCA(j),j);                          % [-] (M,1) Expanded part of the miss distance
    mdConst(j) = norm(posConst(:,pp.NCA(j)) - ...
                                pp.secondary(j).cart(1:3,pp.NCA(j)))^2;    % [-] (M,1) Constant part of the miss distance
end
valRelPos = reshape(mdMaps',3*M,1);                                        % [-] (3M,1) Values of the coefficients of the relative position
valVc     = ones(M,1);                                                     % [-] (M,1)   Values of the coefficients of the virtual buffer

A = sparse([indRelPos1; indVc1], [indRelPos2; indVc2], [valRelPos; valVc]);% [-] (M,*) Build the sparse matrix                    
A = [A sparse(1,nOpt-size(A,2))];                                          % [-] (6N,nOpt) Add zeros to have the correct number of columns

b_lo    = pp.lim + md_exp - mdConst;                                       % [-] (M,1) Residual, comprehensive of the rDeb term to subtract the position of the secondary
b_up    = inf(M,1);                                                        % [-] (M,1) Lower limit of the constraint
end