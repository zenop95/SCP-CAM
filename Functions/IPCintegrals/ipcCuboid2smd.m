function smdLim = ipcCuboid2smd(P, r, R, Pc_target, drMax, fCub)
% invertIPoC
% -------------------------------------------------------------------------
% Converts a IPC limit into SMD using the local inversion of the cuboid
%
%   minimize    || rhat - r ||
%   subject to  Pc(rhat) = Pc_target
%               r - drMax <= rhat <= r + drMax
%
% where Pc() is the cuboid iPoC function.
%
% INPUTS:
%   r          (3x1) original relative position
%   Pc_target  (scalar) desired cuboid PoC value (P_bar_IC in the paper)
%   P, R       (*) additional parameters required by the cuboid function
%   drMax      (3x1) max deviation for constraints
%   fCub       function handle: Pc = fCub(P, rhat, R)
%
% OUTPUTS:
%   rhat_opt   (3x1) optimized point closest to r satisfying Pc = Pc_target
%   fval       optimal cost ||rhat - r||
%   exitflag   fmincon exit condition
%   output     fmincon output struct
%
% -------------------------------------------------------------------------

%% Bounds: r - drMax <= rhat <= r + drMax
lb = r - drMax;
ub = r + drMax;

%% Initial guess (recommended in paper)
r0 = r;

%% Cost function (Eq. 17a)
obj = @(rhat) norm(rhat - r);

%% Nonlinear equality constraint
nonlcon = @(rhat) deal([], cuboidIpc(rhat,P,R) - Pc_target);

%% Optimization options
opts = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'off', ...
    'SpecifyObjectiveGradient', false, ...
    'SpecifyConstraintGradient', false, ...
    'MaxFunctionEvaluations', 2000);

%% Solve
r_opt = fmincon(obj, r0, [], [], [], [], lb, ub, nonlcon, opts);
smd = r'*(P\r);
smdLim = r_opt'*(P\r_opt);
end