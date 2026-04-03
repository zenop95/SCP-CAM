function [A,b_lo,b_up] = homConstr(N,m,sl,uP,k,G_db,uMin,uMax,pp)
% smdGrdContr computes linear matrix and the constant term for the SMD
% gradient constraint in MOSEK
%
% INPUT: N    = [-] Number of nodes in the problem. 
%        m    = [-] Number of optimization variables per node
%        xVec = [-] m*N previous iteration optimization vector
%        k    = [-] Homotopy constant
%        uMin = [-] Minimum constrol acceleration
%        uMax = [-] Maximum constrol acceleration
% 
% OUTPUT: A     = Linear matrix of coefficients
%         b     = Known terms vector
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------   
nOpt  = m*N+sl;    
A = zeros(2*N,nOpt);
b_up = zeros(2*N,1);
b_lo = zeros(2*N,1);
for i = 1:N
   A(i,pp.index.ctrlCone + m*(i-1)) = 1;
   A(i,pp.index.homotopy + m*(i-1)) = -Rmib(uP(i),k,uMin,uMax);
end
% for i = 1:N
%    A(i+N,pp.index.homotopy + m*(i-1)) = dRmib(uP(i),k,uMin,uMax);
%    b_up(i+N) = G_db - Rmib(uP(i),k,uMin,uMax);
%    b_lo(i+N) = -inf;
% end
A = sparse(A);
end


function res = Rmib(u,k,uMin,uMax)
    gMax = uMax - uMin;
    res  = 1-1/(1+exp(k*(u-uMin)/gMax)) + 1/(1+exp(k*gMax));
end

function res = dRmib(u,k,uMin,uMax)
    gMax = uMax - uMin;
    res  = -k/gMax*exp(k*(u-uMin)/gMax)/(exp(k*(u-uMin)/gMax)+1)^2;
end