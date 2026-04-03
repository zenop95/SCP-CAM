function [pointEci, gradEci, ellipsoids] = pointOnEllipsoid(Dr,P,smdLim,pp)
% pointOnEllipsoid solves the convex sub-problem to find the point on the
% ellipsoid closest to Dr.
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

[semiaxes,cov2eci] = defineEllipsoid(P,smdLim);
DrCov              = (cov2eci'*Dr)./semiaxes;

%% Quadratic optimization problem
% Spherical covariance matrix
% Pcov = scCov*cov2eci'*P*cov2eci*scCov;
prob = struct();
% Specify the linear objective terms.
prob.c = -2*DrCov';
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
% prob.qcval  = [1/Pcov(1,1), 1/Pcov(2,2), 1/Pcov(3,3)];
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
%% Find the normal vector to the tangent plane on the point on the ellipsoid
pointCov = pointCovSc.*semiaxes;
gradCov  = 2*pointCov./semiaxes.^2;
gradEci  = cov2eci*gradCov/norm(gradCov);
pointEci = cov2eci*pointCov;

%% Postprocess output
ellipsoids = struct('semiaxes', semiaxes*pp.scaling(1)/sqrt(smdLim), ...   % [m] (3x1) semiaxes of the ellipsoid
                   'Dr',        Dr*pp.scaling(1), ...                      % [m] (3x1) Relative position of the s/c in ECI
                   'grad',      gradEci, ...                               % [-] (3x1) SMD gradient direction
                   'z',         pointEci*pp.scaling(1), ...                % [m] (3x1) position of the point on the ellipsoid's surface
                   'cov2eci',   cov2eci ...                                % [-] (3x3) DCM from the covariance frame to ECI
                    );
end