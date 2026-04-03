function [] = showEllipsoidsEci(pp)
%showEllipsoidsEci plots the relative trajectories with ellipsoid
%visualization
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
lastIter   = pp.majorIter(end);
M          = length(pp.secondary);
n_conj     = M/pp.gmmOrder;
if pp.gmmOrder > 1
    absMan = nan(3,M);
    absNoMan = nan(3,M);
    for j = 1:M
        absMan(:,j)   = pp.validationAbsTraj(1:3,pp.NCA(j));
        absNoMan(:,j) = pp.ballisticTraj(1:3,pp.NCA(j));
    end
end
for ii = 1:n_conj
    figure()
    for jj = 1:pp.gmmOrder
        j          = jj+(ii-1)*pp.gmmOrder;
        i          = pp.NCA(j);
        P          = lastIter.P;
        smdLim     = lastIter.sqrMahaLim(j);
        pNewV      = absMan(:,j)*pp.scaling(1);
        pOld       = absNoMan(:,j)*pp.scaling(1);
        PEci       = P(:,:,i,j);

        [semiaxes,cov2eci] = defineEllipsoid(PEci,smdLim);
        a          = semiaxes(1)*pp.scaling(1);
        b          = semiaxes(2)*pp.scaling(1);
        c          = semiaxes(3)*pp.scaling(1);
        [ax,ang]   = quat2axang(dcm2quat(cov2eci));
        center     = pp.secondary(j).cart(1:3,i)*pp.scaling(1);
        [X,Y,Z]    = ellipsoid(center(1),center(2),center(3),a,b,c,100);
        S          = surf(X,Y,Z,'FaceAlpha',0.3);
        shading interp
        rotate(S,ax,-ang*180/pi,center)
        hold on
        if jj == 1
            plot3(pOld(1),pOld(2),pOld(3),'o','LineWidth',2)
            plot3(pNewV(1),pNewV(2),pNewV(3),'o','LineWidth',2)
        end
    end
    legend('d_m^2 Ellipse','Original','Optimised')
    xlabel('x_{ECI} [km]'); ylabel('y_{ECI} [km]'); zlabel('z_{ECI} [km]')
    grid on
    hold off
end
end