function [] = showEllipseBplane(pp)
%relTrajEllipsoids plots the relative trajectories with ellipsoid
%visualization
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
lastIter = pp.majorIter(end);
j          = 1;
i          = pp.NCA(j);
P          = lastIter.P;
smdLim     = lastIter.sqrMahaLim(j);
smdLimOld  = pp.majorIter(2).sqrMahaLim(j);
pNew       = pp.validationTraj(1:3,i,j)*pp.scaling(1);
pOld       = pp.ballisticRelTraj(1:3,i,j)*pp.scaling(1);
PEci       = P(:,:,i,j);
e2b        = pp.e2b([1 3],:,j);
PB         = e2b*PEci*e2b';

[semiaxes,cov2b] = defineEllipsoid(PB,smdLim);
a          = semiaxes(1)*pp.scaling(1);
b          = semiaxes(2)*pp.scaling(1);
t          = 0:0.001:2*pi;
x          = a*cos(t);
y          = b*sin(t);
ellCov     = [x; y];
ellB       = nan(2,length(t));
off = zeros(2,1);
%         off = e2b*pOld;
for k = 1:length(t)
    ellB(:,k) = cov2b*ellCov(:,k) - off;
end
[semiaxes,cov2b] = defineEllipsoid(PB,smdLimOld);
a          = semiaxes(1)*pp.scaling(1);
b          = semiaxes(2)*pp.scaling(1);
t          = 0:0.001:2*pi;
x          = a*cos(t);
y          = b*sin(t);
ellCov     = [x; y];
ellBOld    = nan(2,length(t));
for k = 1:length(t)
    ellBOld(:,k) = cov2b*ellCov(:,k) - off;
end
pNewB    = e2b*pNew - off;
pOldB    = e2b*pOld - off;
figure()
plot(ellB(1,:),ellB(2,:),'k','HandleVisibility','off');
hold on    
plot(pOldB(1),pOldB(2),'k','Marker','diamond','HandleVisibility','off')
plot(pNewB(1),pNewB(2),'o','Color',[0 0.4470 0.7410],'LineWidth',2)
grid on
box on
xlabel('$\xi$ [km]'); ylabel('$\zeta$ [km]')
axis equal
hold off
end