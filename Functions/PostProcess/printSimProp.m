function [] = printSimProp(pp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
disp('Simulation run paramters');
disp(['Orbit: ', pp.orbit]);
disp(['Initial ephemeris time: ', num2str(pp.et), ' s']);
disp(['Simulation final time: ', num2str(pp.tf/pp.T), ' orbits']);
disp(['Simulation time step: ', num2str(pp.dt), ' s']);
disp(['Gravity harmonics: ', num2str(pp.aida.gravOrd)]);
disp(['Atmosphere: ', num2str(pp.aida.flag1)]);
disp(['SRP: ', num2str(pp.aida.flag2)]);
disp(['Third body: ', num2str(pp.aida.flag3)]);
disp(['Maximum thrust: ', num2str(pp.maxThrust), ' N']);
disp(['CA metric: ', pp.obj]);
disp(['CA metric limit: ', num2str(pp.lim)]);
disp(['Initial CA time: ', num2str(pp.NCA0*pp.dt/pp.T), ' orbits']);
disp(['Final CA time: ', num2str(pp.NCAf*pp.dt/pp.T), ' orbits']);
disp(['Minimum \Delta V constraint: ', num2str(pp.enableHomotopy)]);
disp(['Station keeping: ', num2str(pp.stationKeeping)]);
disp(['Station keeping target: ', num2str(pp.enableSkTarget)]);
if pp.stationKeeping
    disp(['Initial SK time: ', num2str(pp.NSK0*pp.dt/pp.T), ' orbits']);
    disp(['Final SK time: ', num2str(pp.NSKf*pp.dt/pp.T), ' orbits']);
end
disp(['SMD sensitivity constraint: ', num2str(pp.enableSmdGradConstraint)]);
if pp.enableSmdGradConstraint
    disp(['Max IPC deviation at HBR distance: ', num2str(pp.maxIpcDeviation*pp.lim)]);
end
if pp.enableSkTarget
    disp(['Target weight in objective function: ', num2str(pp.targW)]);
end
end