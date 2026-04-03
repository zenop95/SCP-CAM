function [] = showCylinders(pp)
%relTrajEllipsoids plots the relative trajectories with ellipsoid
%visualization
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
colors = [0 0.4470 0.7410;
                  0.8500 0.3250 0.0980;
                  0.9290 0.6940 0.1250;
                  0.4940 0.1840 0.5560;
                  0.4660 0.6740 0.1880;
                  0.3010 0.7450 0.9330;
                  0.6350 0.0780 0.1840];
lastIter   = pp.majorIter(end);
M          = length(pp.secondary);
n_conj     = pp.n_conj;
if pp.gmmOrder > 1
    absMan = nan(3,M);
    absNoMan = nan(3,M);
    for j = 1:M
        absMan(:,j)   = pp.validationAbsTraj(1:3,pp.NCA(j));
        absNoMan(:,j) = pp.ballisticTraj(1:3,pp.NCA(j));
    end
end
faceAlpha  = log10(pp.lim./pp.lims); faceAlpha = faceAlpha/max(faceAlpha);
faceAlpha = 0*ones(M,1);
for ii = 1:n_conj
    figure()
    h = legend('show','location','best');
    for jj = 1:pp.gmmOrder
        j          = jj+(ii-1)*pp.gmmOrder;
        i          = pp.NCA(j);
        P          = lastIter.P;
        smdLim     = lastIter.sqrMahaLim(j);
        pNewV      = absMan(:,j)*pp.scaling(1);
        pOld       = absNoMan(:,j)*pp.scaling(1);
        PEci       = P(:,:,i,j);
        e2b        = pp.e2b(:,:,j);
        e2b2       = e2b;
        e2b2(2,:)  = [];
        PB2d       = e2b2*PEci*e2b2';
    
        [semiaxes1,cov2b2d] = defineEllipsoid(PB2d,smdLim);        
        a          = semiaxes1(1)*pp.scaling(1);
        b          = semiaxes1(2)*pp.scaling(1);
        L          = 50;
        n          = 40;
        center     = pp.secondary(j).cart(1:3,i)*pp.scaling(1);
        v          = linspace(-L/2,L/2,2); 
        t          = linspace(0,2*pi,n);
        [U,V]      = meshgrid(t,v);
        X          = a*cos(U) + center(1); 
        Y          = b*sin(U) + center(2);
        Z          = V        + center(3);
        S          = surf(X,Y,Z,'FaceAlpha',faceAlpha(j),'HandleVisibility','off');
%         shading interp
        cov2e = e2b2'*cov2b2d;
        cov2e3d = cov2e;
        cov2e3d(:,3) = [sqrt(1-norm(cov2e(1,:))^2);
                        sqrt(1-norm(cov2e(2,:))^2)
                        sqrt(1-norm(cov2e(3,:))^2)]';
        if det(cov2e3d) < 0; cov2e3d(:,3) = -cov2e3d(:,3); end
        [ax,ang]   = quat2axang(dcm2quat(cov2e3d));
        rotate(S,ax,-rad2deg(ang),center)
        hold on
        plot3(pOld(1),pOld(2),pOld(3),'diamond','LineWidth',2,'color',colors(jj,:), ...
            'DisplayName',['Original', num2str(jj)])
        plot3(pNewV(1),pNewV(2),pNewV(3),'o','LineWidth',2,'color',colors(jj,:), ...
            'DisplayName',['Optimized', num2str(jj)])
        ellCov = [a*cos(t);b*sin(t)];
        ellEci = nan(3,n);
        for k = 1:n
            ellEci(:,k) = cov2e*ellCov(:,k) + center;
        end
        hold on
        name = strcat('B-plane ',num2str(jj));
        plot3(ellEci(1,:),ellEci(2,:),ellEci(3,:), ...
                '-.','color',colors(jj,:),'LineWidth',2,'DisplayName',name)
    end
    legend
    xlabel('$x_{ECI}$ [km]'); ylabel('$y_{ECI}$ [km]'); zlabel('$z_{ECI}$ [km]')
    grid on
    hold off
end
end