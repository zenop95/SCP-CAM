function prob = createProblem(m,N,s,pp)
% linConvexProblem Solves the linear convex optimization problem of the
% multinode trajectory optimization using the conic optimization tool of 
% MOSEK.
%
% INPUT: m    = [-] (1,1)  Number of optimization variables per node
%        N    = [-] (1,1) Number of nodes in the problem. 
%        s    = [struct] Structure containing the SCVX matrices and vectors
%        pp   = [struct] Postprocess structure
% 
% OUTPUT: prob: [struct] Structure of the problem as requested by MOSEK
%
% Documentation: https://docs.mosek.com/9.2/toolbox/tutorial-cqo-shared.html 
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
s.b_gradBound_lo = zeros(length(s.b_gradBound),1);

%% Specify limits and cost coefficients
prob.blx = s.limsLo;
prob.bux = s.limsUp;
prob.c   = s.c;

%% Specify constraint coefficient matrix
prob.a   = [s.A_dx; 
            s.A_dyn; 
            s.A_rel; 
            s.A_ca; 
            s.A_sk;
            s.A_alt
            s.A_skTar; 
            s.A_grad;
            s.A_gradBound;
            s.A_hom];

%% Specify constraint upper and lower bounds
prob.blc = [s.b_dx_lo; 
            s.b_dyn; 
            s.b_rel;    
            s.b_ca_lo; 
            s.b_sk_lo; 
            s.b_alt_lo; 
            s.b_skTar; 
            s.b_grad;  
            s.b_gradBound_lo
            s.b_hom_lo];

prob.buc = [s.b_dx_up; 
            s.b_dyn;  
            s.b_rel;    
            s.b_ca_up; 
            s.b_sk_up;
            s.b_alt_up; 
            s.b_skTar; 
            s.b_grad;
            s.b_gradBound;
            s.b_hom_up];

%% Specify the cones
prob.cones = buildCones(N,m,pp);
end