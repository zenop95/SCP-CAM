function cones = buildCones(N,m,pp)
% BuildCones Builds the cones as required by the conic optimization tool of 
% MOSEK.
%
% INPUT: N    = [-] (1,1) Number of nodes in the problem. 
%        m    = [-] (1,1)  Number of optimization variables per node
%        pp   = [struct] Postprocess structure
% 
% OUTPUT: cones: [struct] Structure of the cones as requested by MOSEK
%
% Documentation: https://docs.mosek.com/9.2/toolbox/tutorial-cqo-shared.html 
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

%% Initialize variables
i        = pp.index;
ind      = reshape(1:m*N,m,N);
typeGrad = []; conesGrad = []; conesSubptrGrad = [];
typeVc   = []; conesVc   = []; conesSubptrVc   = [];

%% Acceleration cones
l             = length(i.ctrl) + 1;
typeAcc       = zeros(1,N);
conesAcc      = [reshape(ind([i.ctrlCone,i.ctrl],:),1,l*N)];
conesSubptrAcc = 1:l:l*N; % starting indeces of cone

%% SMD gradient cones
if pp.enableSmdGradConstraint
    typeGrad        = zeros(1,N);
    conesGrad       = reshape(ind([i.gradCone,i.grad],:),1,4*N);
    conesSubptrGrad = l*N+1:4:(4+l)*N; % starting indeces of cones
end

%% Virtual control
off           = pp.enableSmdGradConstraint*4+l;
typeVc        = zeros(1,N);
indVc         = ind([i.vcCone,i.vc],:);
s             = size(indVc,1);
conesVc       = reshape(indVc,1,s*N);
conesSubptrVc = off*N+1:s:(off + s)*N;


%% Build cones structure
cones.type   = [typeAcc,        typeGrad,        typeVc];
cones.sub    = [conesAcc,       conesGrad,       conesVc];
cones.subptr = [conesSubptrAcc, conesSubptrGrad, conesSubptrVc];

end


