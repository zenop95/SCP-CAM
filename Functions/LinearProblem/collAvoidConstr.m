function [A,b_up,b_lo,minorIter] = collAvoidConstr(m,nOpt,N0,Nf,rel_pos, ...
                                              P,smdLim,e2b,j,minorIter,pp)
% CollAvoidConstr computes linear matrix and the constant term for the 
% collision avoidance constraint in terms of linearized keep-out zone in 
% MOSEK. The nodewise vector constraint implemented is the following
% 
%               Nabla(dm^2)(z)*(dr-z) >= 0
%
% INPUT: m       = [-] (1,1)   Number of optimization variables per node
%        nOpt    = [-] (1,1)   Total number of optimization variables. 
%        N0      = [-] (1,1)   Initial node for the collision avoidance. 
%        Nf      = [-] (1,1)   Final node for the collision avoidance. 
%        rel_pos = [-] (6,N)   Previous iteration vector of expansion
%                              points, with the points inside the ellipsoid   
%                              shifted to the surface of the ellipsoids
%        P       = [-] (3,3,N) Covariance matrices
%        smdLim  = [-] (N,1)   Vector of SMD limit
%        pp      = [struct]    Postprocess structure
% 
% OUTPUT: A         = [-] (Nf-N0+1,nOpt) Linear matrix of coefficients
%         b_up      = [-] (Nf-N0+1,1)    Upper limit for the constraint
%         b_lo      = [-] (Nf-N0+1,1)    Lower limit for the constraint
%         minorIter = [struct]           Structure of minor iteration
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------   
if pp.enableBplaneAvoidance
    b_up = inf;
    A    = zeros(1,nOpt);
    i    = pp.NCA(j);
    [z, grad, ellipsoids] = pointOnEllipseBplane(rel_pos(:,i), ...
                                             P(:,:,i),smdLim,e2b,pp);      % Solve the proejction problem to find he closest point on the surface of the ellipsoid
    b_lo = dot(grad,z + pp.secondary(j).cart(1:3,i));                      % Define the linearization of the CA constraint
    A(pp.index.state(1:3) + m*(i-1)) = grad';
    A(pp.index.caVc(j) + m*(i-1)) = 1;                                     % Define the virtual buffer coefficients
    minorIter.gradCA(:,i) = grad;
    minorIter.z(:,i)      = z;
    minorIter.ellipsoid(i) = ellipsoids;
else
    b_up = inf(Nf-N0+1,1);
    b_lo = zeros(Nf-N0+1,1);
    A    = zeros(Nf-N0+1,nOpt);
    for i = N0:Nf
        [z, grad, ellipsoids] = pointOnEllipsoid(rel_pos(:,i), ...
                                            P(:,:,i),smdLim(i),pp);        % Solve the proejction problem to find he closest point on the surface of the ellipse
        b_lo(i+1-N0) = dot(grad,z + pp.secondary(j).cart(1:3,i));          % Define the linearization of the CA constraint
        A(i+1-N0, pp.index.state(1:3) + m*(i-1)) = grad';
        A(i+1-N0, pp.index.caVc(j) + m*(i-1)) = 1;                         % Define the virtual buffer coefficients
        minorIter.ellipsoid(i) = ellipsoids;
        minorIter.gradCA(:,i) = grad;
        minorIter.z(:,i)      = z;
    end    
end
A = sparse(A);
end