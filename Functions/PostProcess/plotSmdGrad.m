function plotSmdGrad(t,pp)
%PlotSmdGrad plots the original and optimised evolution of the gradient of
% the SMD
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
N = pp.N;
gradOld   = pp.ballisticSmdGrad;
gradNew   = pp.smdGrad;
normOld = normOfVec(gradOld)/pp.scaling(1);
normNew = normOfVec(gradNew)/pp.scaling(1);
figure()
semilogy(t,normOld)
hold on
semilogy(t,normNew)
try 
    plot(t,pp.majorIter(end).gradLim/pp.scaling(1),'k--')
catch 
end
xlabel('Number of orbits [-]')
ylabel('$|\nabla d_M^2| [\mathrm{m}^{-1}]$')
axis tight
ylim([min(min(normNew),min(normOld))/10, max(max(normNew),max(normOld))*10])
legend('Original','Optimised')
grid on
hold off

figure()
semilogy(t,pp.deltaIpcLin)
hold on
semilogy(t,pp.deltaIpc)
semilogy(t,pp.ipcLim*pp.maxIpcDeviation*ones(N,1),'k--')
axis tight
ylim([1e-12 pp.ipcLim*100])
title('P_{IC} difference at HBR distance in \nabla d_m^2 direction')
grid on
xlabel('Number of orbits [-]')
ylabel('$\Delta P_{IC}$ [-]')
legend('Shifted computed at HBR distance','Approximated HBR distance', ...
    'Limit ($\Delta P_{IC}$)','Interpreter','Latex')

% figure()
% semilogy(t,pp.ipcAtOffset)
% hold on
% semilogy(t,pp.ipcLin)
% semilogy([t(1),t(end)],pp.ipcLim*pp.lim*ones(2,1),'k--')
% axis tight
% ylim([1e-10 pp.ipcLim*100])
% title('P_{IC} at HBR distance in \nabla d_m^2 direction')
% grid on
% xlabel('t/T [-]')
% ylabel('P_{IC} [-]')
% legend('Shifted computed at HBR distance','Approximated HBR distance', ...
%     'Limit (\Delta P_{IC})')

end