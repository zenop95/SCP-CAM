function [] = dispResults(pp)
%DispResults display some accuracy results of the simulation
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
nMinor = 0;
for i = 1:length(pp.majorIter)
    nMinor = nMinor + length(pp.majorIter(i).minorIter);
end
disp(['The SCVX algorithm took: ', num2str(pp.simTime-pp.timeSubtr), ' s']);
disp(['Total DV: ', num2str(pp.DvTot), ' m/s']);
disp(['Last major Iteration error: ', num2str(pp.majorIter(end).err)]);
disp(['Number of major iterations: ', num2str(length(pp.majorIter))]);
disp(['Number of minor iterations: ', num2str(nMinor)]);
disp(['Maximum validation error: ', num2str(pp.valErr), ' km']);
if pp.enableSkTarget
    disp(['Error in final target position: ', ...
        num2str(norm(pp.validationAbsTraj(1:3,end) - ...
                     pp.cartTarget(1:3))*pp.scaling(1)),' km']);
    disp(['Error in final target velocity: ', ...
        num2str(norm(pp.validationAbsTraj(4:6,end) - ...
                     pp.cartTarget(4:6))*pp.scaling(4)),' km/s']);
end
if pp.fastEncounter
   disp(['PoC type for optimization: ', pp.PoCType])
   disp(['Maximum PoC: ', num2str(pp.PoCTotMaximum)])
   disp(['Alfriend PoC: ', num2str(pp.PoCTotConstant)])
   disp(['Chan PoC: ', num2str(pp.PoCTotChan)])
   disp(['Alfano PoC: ', num2str(pp.PoCTotAlfano)])
else
   disp(['IPoC type for optimization: ', pp.ipc_type])
   disp(['Maximum value of IPoC: ', num2str(max(pp.ipcMan))])
end
end