function [] = plotSingleEllipsoid(majorIter,node)
%PlotEllipsoidPoints plots the covariance ellipsoid for the required node
%and the interior point or the successive convexification points.
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
lll = [];
leg = 'ell';
% if isfield("ellipsoid",majorIter.minorIter(1))
    cov2eci  = majorIter.minorIter(1).ellipsoid(node).cov2eci;
    a       = majorIter.minorIter(1).ellipsoid(node).semiaxes(1);
    b       = majorIter.minorIter(1).ellipsoid(node).semiaxes(2);
    c       = majorIter.minorIter(1).ellipsoid(node).semiaxes(3);
    [ax,ang] = quat2axang(dcm2quat(cov2eci));
    [X,Y,Z] = ellipsoid(0,0,0,a,b,c,1000);
    figure()
    S = surf(X,Y,Z,'FaceAlpha',0.3);
    shading interp
    rotate(S,ax,-ang*180/pi)
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    hold on    
    for j = 1:length(majorIter.minorIter)
        point   = majorIter.minorIter(j).ellipsoid(node).z;
        Dr      = majorIter.minorIter(j).ellipsoid(node).Dr;
        plot3([point(1) Dr(1)],[point(2) Dr(2)],[point(3) Dr(3)],'o')
        axis equal
        lll = [lll, majorIter.minorIter(j).ellipsoid(node).Dr, ...
                    majorIter.minorIter(j).ellipsoid(node).z];
        leg = [leg; string(['Iteration', num2str(j)])];
    end
%     xu = max(lll(1,:));    yu = max(lll(2,:));    zu = max(lll(3,:));
%     xl = min(lll(1,:));    yl = min(lll(2,:));    zl = min(lll(3,:));
%     xlim([xl-.1,xu+.1]); ylim([yl-.1,yu+.1]); zlim([zl-.1,zu+.1])
%     legend(leg)
    hold off
% end
