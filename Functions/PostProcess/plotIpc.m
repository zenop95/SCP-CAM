function plotIpc(t, pp)
%PlotIpc Plots the orginal and optimised IPC evolutions
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
figure()
semilogy(t,pp.ballisticIpcTot)
hold on
semilogy(t,pp.ipcMan)
semilogy([t(1),t(end)],pp.ipcLim*ones(1,2),'k--')
xlabel('Number of orbits [-]')
ylabel('$P_{TIC}$ [-]')
legend('Ballistic','Optimized','P_{IC} limit')
axis tight
ylim([pp.ipcLim/1000,pp.ipcLim*1000])
grid on
hold off
end