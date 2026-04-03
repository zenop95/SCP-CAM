function [pointEci, gradEci, ellipsoids] = ...
                              pointOnEllipseBplane(Dr,PEci,smdLim,e2b,pp)
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
[cov2eci,~] = eig(PEci);
if det(cov2eci) < 0; cov2eci(:,1) = -cov2eci(:,1); end
%pass in 2d to solve B-plane problem
e2b2      = e2b;
e2b2(2,:) = [];
rB        = e2b2*Dr;
PB        = e2b2*PEci*e2b2';

[semiaxes,cov2b] = defineEllipsoid(PB,smdLim);
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
prob.qosubi = [1 2]';
prob.qosubj = [1 2]';
prob.qoval  = [2 2]';

% Specify the quadratic terms of the constraints.
% Here there is also the index k because it represents the constraint
% number, if you have more than one constraint you start adding twos.
prob.qcsubk = [1 1];
prob.qcsubi = [1 2];
prob.qcsubj = [1 2];
prob.qcval  = ones(1,2);

% Specify the linear constraint matrix (no constraint)
prob.a      = zeros(1,2);

% Specify the bounds of the constraint
prob.buc    = 1/2;

% Specify the bounds of the decision variables
prob.bux    =  ones(2,1);
prob.blx    = -ones(2,1);

% Solve the problem
[~,res]    = mosekopt('minimize echo(0)',prob);
pointCovSc = res.sol.itr.xx; % expressed in cov
pointCov   = pointCovSc.*semiaxes;

%% Find the normal vector to the tangent plane on the point on the ellipsoid
gradCov = 2*pointCov./semiaxes.^2;
gradB   = cov2b*gradCov; 
pointB  = cov2b*pointCov; 

% pass in 3D again
gradEci  = e2b2'*gradB/norm(gradCov);
pointEci = e2b2'*pointB;
%% Postprocess output
ellipsoids = struct('semiaxes', semiaxes*pp.scaling(1)/sqrt(smdLim), ...   % [m] (3x1) semiaxes of the ellipsoid
                   'Dr',        rB*pp.scaling(1), ...                      % [m] (3x1) Relative position of the s/c in ECI
                   'grad',      gradEci, ...                               % [-] (3x1) SMD gradient direction
                   'z',         pointEci*pp.scaling(1), ...                % [m] (3x1) position of the point on the ellipsoid's surface
                   'cov2Eci',   cov2eci, ...                               % [-] (3x3) DCM from the covariance frame to ECI
                   'cov2b',     cov2b ...                                  % [-] (3x3) DCM from the covariance frame to the B-plane
                    );
end