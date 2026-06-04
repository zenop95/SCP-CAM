%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                                                         %
% Sequential Convex Program for Collision Avoidance Manoeuvres (SCP CAM)  %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%--------------------------------------------------------------------------
%        ███████╗ ██████╗██████╗      ██████╗ █████╗ ███╗   ███╗          L
%        ██╔════╝██╔════╝██╔══██╗    ██╔════╝██╔══██╗████╗ ████║          L
%        ███████╗██║     ██████╔╝    ██║     ███████║██╔████╔██║          L
%        ╚════██║██║     ██╔═══╝     ██║     ██╔══██║██║╚██╔╝██║          L
%        ███████║╚██████╗██║         ╚██████╗██║  ██║██║ ╚═╝ ██║          L
%         ══════╝ ╚═════╝╚═╝          ╚═════╝╚═╝  ╚═╝╚═╝     ╚═╝          L
%--------------------------------------------------------------------------
%
% This project implements Sequential Convex Programming (SCP) methods for
% Collision Avoidance Manoeuvres (CAMs). This is a complete software to
% compute CAMs for any kind of conjunctions, both short- and long-term,
% with single or multiple secondaries.
% 
% This script initializes the environment, 
% configures simulation properties, and runs the major-iteration loop 
% for trajectory optimization using successive convexification. 
% It optionally applies linear IPC/PoC constraints, adapts SMD limits, 
% validates the final solution, and runs postprocessing.
%
% High-level Steps
%  1. initializePath - set up MATLAB path and toolboxes required by the project
%  2. Configure simulation properties for the chosen orbit (LEO/ GEO)
%  3. Prepare Gaussian mixture models, firing nodes, and AIDA assets
%  4. Run majorSolve to perform major iterations and update trajectory
%  5. Optionally embed linear PoC constraints and re-run
%  6. Adapt SMD limits iteratively if enabled
%  7. Validate results and run postprocessing/plotting
%
% Required Inputs / Preconditions (fields in pp)
%  - pp.orbit            : 'leo' or 'geo'
%  - pp.N                : number of discretization nodes
%  - pp.iterMaxMaj, tolMaj: major iteration stopping criteria
%  - pp.gmmOrder, primary, secondary, etc. for encounter modeling
%  - Modification of the simProperties function
%  - functions on path: simProperties, splitInitial, defIndConv,
%                       findFiringNodes, majorSolve,
%                       validate_pp (or validate), selectFigs, mainPostProcess
%
% Outputs
%  - pp.majorIter populated with per-major-iteration data
%  - Figures created by postprocessing and files if downstream routines save them
%
% Compatibility:
%  - MATLAB R2022a and newer (tested up to R2025b)
%  - MOSEK 10.0
%
% Author: Zeno Pavanello
% E-mail: zpav176@aucklanduni.ac.nz
% Date:   2022-2024
%--------------------------------------------------------------------------

clc; clear; %close all;

% 0) Path init / toolboxes
initializePath();

% 1) Build full configuration
pp = simProperties(1);

% 2) Run SCP pipeline
pp = runSCP(pp);

% 3) Validation
pp = validateSolution(pp);

% 4) Postprocessing
pp.figs = selectFigs(pp);
mainPostProcess(pp);
