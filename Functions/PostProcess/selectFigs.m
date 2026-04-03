function figs = selectFigs(pp)
% SelectFigs defines which figures to show in the post process
%   Detailed explanation goes here
figs.sqMaha           = false;
figs.ipc              = true;
figs.eci              = false;
figs.ecef             = false; 
figs.lla              = false;
figs.alt              = false;
figs.ellipsoids       = false;
figs.singleSphere     = false;
figs.minorConvergence = false;
figs.majorConvergence = false;      
figs.relTraj          = false;
figs.relTrajDv        = false;
figs.relTrajEll       = false;
figs.relTrajEachMajor = false;
figs.relTrajEachMinor = false;   
figs.linCaConstr      = false;   
figs.dv               = true;
figs.acc              = false;
figs.grad             = false;
figs.gamma            = false;
figs.nli              = false;
figs.limsEvolution    = false;
figs.PoCTime          = false;
if pp.fastEncounter
    figs.ipc          = false;
    figs.smd          = false;
    figs.grad         = false;
end
end