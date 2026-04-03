function [] = showMultipleEllipsoids(pp,i)
%showMultipleEllipsoids plots the relative trajectories with ellipsoid
%visualization
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
nn         = 1000;
lastIter   = pp.majorIter(end);
M          = length(pp.secondary);
color1     = [0.4660 0.6740 0.1880];
color2     = [0 0.4470 0.7410];
color3     = [0.8500 0.3250 0.0980];
limLog     = -log10(pp.lims);
limMin     = 10;
limMax     = 6;
limCode    = linspace(limMin,limMax,nn);
colorCode  = color1 + 2/nn*(color2-color1).*repmat([0:nn/2,nn/2*ones(1,nn/2)],3,1)' + ...
                2/nn*(color3-color2).*repmat([zeros(1,nn/2),0:nn/2],3,1)';
[~,colInd] = min(abs(repmat(limCode,M,1)-limLog),[],2);
colors     = colorCode(colInd,:);

figure()
hold on
rN             = pp.validationAbsTraj(1:3,i)*pp.scaling(1);
rO             = pp.ballisticTraj(1:3,i)*pp.scaling(1);
for j = 1:M
    offCenter      = pp.secondary(j).cart(1:3,i)*pp.scaling(1);
    smdLim         = pp.majorIter(end).sqrMahaLim(i,j);
    PEci           = lastIter.P(:,:,i,j);
    [semiax,cov2e] = defineEllipsoid(PEci,smdLim);
    semiax = semiax*pp.scaling(1);
    [ax,ang] = quat2axang(dcm2quat(cov2e));
    [X,Y,Z] = ellipsoid(offCenter(1),offCenter(2),offCenter(3), ...
                        semiax(1),semiax(2),semiax(3),100);
    S       = surf(X,Y,Z,'FaceAlpha',.1,'LineStyle','none', ...
                'FaceColor',colors(j,:),'HandleVisibility','off');
    shading interp
    rotate(S,ax,-ang*180/pi,offCenter)
    plot3(offCenter(1),offCenter(2),offCenter(3), ...
          'k*','LineWidth',2,'HandleVisibility','off')
end
plot3(rO(1),rO(2),rO(3),'diamond','Color', ...
    [0 0.4470 0.7410],'LineWidth',1,'DisplayName','Ballistic')
plot3(rN(1),rN(2),rN(3),'o','Color', ...
    [0 0.4470 0.7410],'LineWidth',1.5,'DisplayName','Optimized')
xlabel('$x$ [-]')
ylabel('$y$ [-]')
zlabel('$z$ [-]')
grid on
view(3)
legend show
hold off
end