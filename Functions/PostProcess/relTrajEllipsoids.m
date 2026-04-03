function [] = relTrajEllipsoids(pp)
%relTrajEllipsoids plots the relative trajectories with ellipsoid
%visualization
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
firstIter = pp.majorIter(1);
lastIter  = pp.majorIter(end);
N         = pp.N;
pNew      = lastIter.minorIter(end).p;
pOld      = firstIter.pNoMan;
rS        = lastIter.absPMan;
rD        = pp.x_dProp(1:3,:);
ellipsoids = lastIter.minorIter(end).ellipsoid;
smdLim     = lastIter.sqrMahaLim;

figure()
plot3(rD(1,:),rD(2,:),rD(3,:),'.','LineWidth',3)
hold on
plot3(rS(1,:),rS(2,:),rS(3,:),'.','LineWidth',3)
plot3(rS(1,1),rS(2,1),rS(3,1),'o')
plot3(rS(1,38),rS(2,38),rS(3,38),'*')
plot3(rS(1,end),rS(2,end),rS(3,end),'x')
for i = pp.NCA0:1:pp.NCAf
    if smdLim > 1e-4
        axes = [ellipsoids(i).a;ellipsoids(i).b;ellipsoids(i).c]*sqrt(smdLim(i));
        [ellX,ellY,ellZ] = ellipsoid(rD(1,i),rD(2,i),rD(3,i), ...
                                     axes(1),axes(2),axes(3),10);
        S = surf(ellX,ellY,ellZ);%,'EdgeColor','none');
        [ax,ang] = quat2axang(dcm2quat(ellipsoids(i).cov2CW));
        rotate(S,ax,-ang*180/pi,rD(:,i))
    end
end
legend('Secondary','Primary')
xlabel('[m]')
ylabel('[m]')
zlabel('[m]')
% zlim([-300,100])
grid on
hold off

% figure()
% xlabel('[m]')
% ylabel('[m]')
% zlabel('[m]')
% for i = pp.NCA0:pp.NCAf
% 
%     plot3(pOld(1,i),pOld(2,i),pOld(3,i),'-o','LineWidth',2,'color',[0 0.4470 0.7410])
%     hold on
%     title(['Relative trajectory at t/T = ', num2str(pp.t(i)/pp.T)])
%     plot3(pNew(1,i),pNew(2,i),pNew(3,i),'-o','LineWidth',2,'color',[0.8500 0.3250 0.0980])
%         axes = [ellipsoids(i).a;ellipsoids(i).b;ellipsoids(i).c]*sqrt(smdLim(i));
%         [ellX,ellY,ellZ] = ellipsoid(0,0,0,axes(1),axes(2),axes(3),30);
%         [ax,ang] = quat2axang(dcm2quat(ellipsoids(i).cov2CW));
%         S = surf(ellX,ellY,ellZ,'FaceAlpha',1,'LineStyle','none', ...
%                     'FaceColor',[0.9290 0.6940 0.1250]);
%         rotate(S,ax,-ang*180/pi,[0,0,0])
% %         shading interp
%     legend('Original','Optimised','d_m^2 Ellipsoids')
% %     pause
%     alpha(0.1); % increase transparency of objets of previous iterations
% end
% xlabel('x_{ECI}'); ylabel('y_{ECI}'); zlabel('z_{ECI}')
% hold off

% for i = 1:N
%     title('Relative trajectory')
%     axes = [ellipsoids(i).a;ellipsoids(i).b;ellipsoids(i).c]*Lsc;
%     scale = max(axes);
%     dcm  = ellipsoids(i).cov2CW;
% %     rotAxes = dcm*axes;
%     p(:,i) = dcm'*pOld(:,i);
%     r(:,i) = dcm'*pNew(:,i);
% end
% figure() % keep ellipsoid constant and rotate and move the relative position
% xlabel('[m]')
% ylabel('[m]')
% zlabel('[m]')
% plot3(p(1,:),p(2,:),p(3,:),'color',[0 0.4470 0.7410])
% hold on
% plot3(r(1,:),r(2,:),r(3,:),'color',[0.8500 0.3250 0.0980])
% pause
% [ellX,ellY,ellZ] = ellipsoid(0,0,0,axes(1)/scale, ...
%                              axes(2)/scale,axes(3)/scale,30);
% surf(ellX,ellY,ellZ);
% legend('Original','Optimised','Normalized Ellipsoid')
% shading interp
% hold off
% end