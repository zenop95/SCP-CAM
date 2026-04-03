function [] = showSingleCircle(pp)
%relTrajEllipsoids plots the relative trajectories with ellipsoid
%visualization
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
lastIter = pp.majorIter(end);
M        = length(pp.NCA);
lims     = pp.lims;
for v = 1:pp.n_conj
    for jj = 1:pp.gmmOrder
        j          = jj + (v-1)*pp.gmmOrder;
        i          = pp.NCA(j);
        P          = lastIter.P;
        smdLim     = lastIter.sqrMahaLim(j);
        pNew       = pp.validationTraj(1:3,i,j);
        pOld       = pp.ballisticRelTraj(1:3,i,j);
        PEci       = P(:,:,i,j);
        e2b        = pp.e2b(:,:,j);
        e2b(2,:)   = []; 
        PB         = e2b*PEci*e2b';   

        [semiaxes,cov2b] = defineEllipsoid(PB,smdLim);
        e2cov = cov2b'*e2b;
        pNewB(:,j) = e2cov*pNew./semiaxes;
        pOldB(:,j) = e2cov*pOld./semiaxes;
    end
end
lims(ceil(pp.gmmOrder/2)) = sum(lims(1:pp.gmmOrder));
y = [];
for j = 1:pp.gmmOrder
    if pp.secondary(j).cdm && j ~= ceil(pp.gmmOrder/2)
        y = [y,j];
    end
end
pNewB(:,y) = [];
pOldB(:,y) = [];
lims(y)    = [];
t          = 0:0.001:2*pi;
figure()
plot(cos(t),sin(t),'k','LineWidth',2);
hold on    
s1 = scatter(pOldB(1,:),pOldB(2,:),[],log10(lims),'d','filled');
s2 = scatter(pNewB(1,:),pNewB(2,:),[],log10(lims),'filled');
s1.SizeData = 60;
s2.SizeData = 60;
colormap jet;
cb = colorbar;
if pp.gmmOrder > 1
    ylabel(cb,'$\mathrm{log_{10}}(\bar{P}_{C,cs})$','FontSize',16,'Interpreter','latex')
else
    ylabel(cb,'$\mathrm{log_{10}}(\bar{P}_{C,s})$','FontSize',16,'Interpreter','latex')
end
grid on
xlabel('$\tilde{\xi}$ [-]'); ylabel('$\tilde{\zeta}$ [-]')
axis equal
hold off
end