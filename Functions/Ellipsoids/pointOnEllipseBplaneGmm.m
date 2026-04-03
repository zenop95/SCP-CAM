function [pointEci, gradEci] = pointOnEllipseBplaneGmm(rVec,P,smdLim,e2b,conj,pp)
% pointOnEllipse solves the convex sub-problem to find the point on the
% B-plane ellipse closest to Dr.
%
% INPUT: Dr     = [-] (3x1) relative distance between the two bodies
%                 expressed in the frame centered on the center of the 
%                 ellipsoid. If Dr is inside the ellipsoid, the value stays
%                 unchanged.
%        P      = [-] (3x3) covariance matrix of the relative distance
%                 expresse in ECI coordinates
%        smdLim = [-] (1x1) Limit value of the SMD
%
% OUTPUT: point     = [-] (3x1)relative distance that ensures the SMD 
%                     constraint is respected.
%
% Note: Validated using the script testPointOnEllipsoid
% reference Mosek: https://docs.mosek.com/9.2/toolbox/tutorial-qo-shared.html

% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
offset  = nan(3,pp.gmmOrder);
semiaxes = nan(3,pp.gmmOrder);
cov2b    = nan(3,3,pp.gmmOrder);
PB       = nan(3,3,pp.gmmOrder);
center   = pp.secondary(ceil(pp.gmmOrder*(conj+0.5))).cart(1:3,:);
for j = 1:pp.gmmOrder
    jj = j + pp.gmmOrder*(conj-1);
    PEci        = P(:,:,jj);
    %pass in 2d to solve B-plane problem
    e2b       = pp.e2b(:,:,jj);
    Dr        = rVec(:,jj);
    rB        = e2b*Dr;
    PB(:,:,j) = e2b*PEci*e2b';
    [semiaxes,cov2b] = defineEllipsoid(PB(:,:,j),smdLim(jj));
    offset(:,j)     = pp.secondary(jj).cart(1:3,:) - center;
end
DrCov            = cov2b'*rB;
DrCovSc          = DrCov./semiaxes;
%Diagonal covariance matrix
prob = struct();
% Specify the linear objective terms.
prob.c = -2*DrCovSc';
% Specify the quadratic terms of the objective.
% qosubi and qosubj list the indices (first vector row index, second vector
% column index) of the elements of the matrix Q in a sparse format. qoval 
% assigns the values to the Q matrix according to the indexes from before.
prob.qosubi = [1 2 3]';
prob.qosubj = [1 2 3]';
prob.qoval  = [2 2 2]';

% Specify the quadratic terms of the constraints.
% Here there is also the index k because it represents the constraint
% number, if you have more than one constraint you start adding twos.
prob.qcsubk = [1 1 1];
prob.qcsubi = [1 2 3];
prob.qcsubj = [1 2 3];
prob.qcval  = ones(1,3);

% Specify the linear constraint matrix (no constraint)
prob.a      = zeros(1,3);

% Specify the bounds of the constraint
prob.buc    = 1/2;

% Specify the bounds of the decision variables
prob.bux    =  ones(3,1);
prob.blx    = -ones(3,1);

% Solve the problem
[~,res]    = mosekopt('minimize echo(0)',prob);
pointCovSc = res.sol.itr.xx; % expressed in cov
pointCov   = pointCovSc.*semiaxes;

%% Find the normal vector to the tangent plane on the point on the ellipsoid
gradCov = 2*pointCov./semiaxes.^2;
gradB   = cov2b*gradCov; 
pointB  = cov2b*pointCov; 

% pass in 3D again
gradEci  = e2b'*gradB/norm(gradCov);
pointEci = e2b'*pointB;
end