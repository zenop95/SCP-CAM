function plotMd(t, pp)
%PlotIpc Plots the orginal and optimised MD evolutions
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
figure()
plot(t,normOfVec(pp.ballisticRelTraj(1:3,:))*pp.scaling(1))
hold on
plot(t,normOfVec(pp.validationTraj(1:3,:))*pp.scaling(1))
plot([t(1),t(end)],pp.mdLim*ones(2,1),'k--')
xlabel('Number of orbits [-]')
ylabel('$MD$ [km]')
legend('Ballistic','Optimized','MD limit')
axis tight
grid on
hold off
end