function [] = relTrajDv(pp)
%RELTRAJDV Plots the orginal and optimised trajectories with the Dv
%indications
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

N          = pp.N;
pNew       = pp.majorIter(end).minorIter(end).p;
pOld       = pp.majorIter(1).pNoMan;
dv         = pp.dv;
dv         = dv/pp.uMax; % scale the dv vector for better visualization

figure()
plot3(pOld(1,:),pOld(2,:),pOld(3,:))
hold on
plot3(pNew(1,:),pNew(2,:),pNew(3,:))
xlabel('[m]')
ylabel('[m]')
zlabel('[m]')
grid on
for i = 1:N
    if norm(dv(:,i))
        quiver3(pNew(1,i),pNew(2,i),pNew(3,i),dv(1,i),dv(2,i),dv(3,i),'r','AutoScale','on')
    end
end
hold off

end

