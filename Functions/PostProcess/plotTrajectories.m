function plotTrajectories(t,pp)
%PlotTrajectories plots graphs of the trajectories
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
N        = pp.N;
M        = length(pp.secondary);
pNew     = pp.majorIter(end).relP;
for j = 1:M
    pOld(:,:,j)     = pp.ballisticTraj(1:3,:) - pp.secondary(j).cart(1:3,:);
end
val      = pp.validationTraj;
pvalLvlh = pp.pLvlh;
p        = pp.figs;
% z        = pp.majorIter(end).minorIter(end).z;
% grad     = pp.majorIter(end).minorIter(end).gradCA;

%% ECI trajectory
if p.eci
    absOld = pOld(:,:,1) + pp.secondary(1).x(1:3,:);
    plot3(absOld(1,:),absOld(2,:),absOld(3,:))
    hold on
    plot3(pp.validationAbsTraj(1,:),pp.validationAbsTraj(2,:),pp.validationAbsTraj(3,:))
    grid on
    hold off
end
%% Relative trajectory
if p.relTraj
    for j = 1:M
        pO = pOld(:,:,j).*pp.scaling(1);
        pval = pvalLvlh(:,:,j)*pp.scaling(1);
        for i = 1:N
            eci2lvlh    = lvlh2eci(pp.secondary(j).cart(1:3,i), ...
                                   pp.secondary(j).cart(4:6,i))';
            pOLvlh(:,i) = eci2lvlh*pO(:,i);
        end
        figure()
        plot3(pOLvlh(1,:),pOLvlh(2,:),pOLvlh(3,:))
        hold on
        plot3(pval(1,:),pval(2,:),pval(3,:))
    end
        legend('Original','Optimised')
        xlabel('$x_{LVLH}$ [km]')
        ylabel('$y_{LVLH}$ [km]')
        zlabel('$z_{LVLH}$ [km]')
        title('Relative trajectories in debris LVLH frame')
        hold off
        grid on
end

%% Every major iteration relative trajectory 
if p.relTrajEachMajor
    light = [0.3010 0.7450 0.9330];
    dark = [0.6350 0.0780 0.1840];
    ncolors = length(pp.majorIter);
    gradColors = @(n,nn) interp1([1/nn 1],[light;dark],n/nn); 
    figure()
    leg = [];
    colororder([jet(44)])
    for i = 1:length(pp.majorIter)
        for j = 1:N
            eci2lvlh    = lvlh2eci(pp.x_dProp(1:3,j), ...
                                   pp.x_dProp(4:6,j))';
            pMan(:,j) = eci2lvlh*pp.majorIter(i).pMan(:,j);
        end
        plot3(pMan(1,:),pMan(2,:),pMan(3,:),'color',gradColors(i,ncolors))
        hold on
        leg = char([leg; string(['iter = ',num2str(i)])]);
    end
    legend(leg);
    xlabel('[m]')
    ylabel('[m]')
    zlabel('[m]')
    hold off
    grid on
end

%% Every major iteration relative trajectory 
if p.relTrajEachMinor
    light = [0.3010 0.7450 0.9330];
    dark = [0.6350 0.0780 0.1840];
    ncolors = 0;
    for i = 1:length(pp.majorIter)
       ncolors = ncolors + length(pp.majorIter(i).minorIter);
    end
    gradColors = @(n,nn) interp1([1/nn 1],[light;dark],n/nn); 
    figure()
    leg = [];
    l   = 0;
    for i = 1:length(pp.majorIter)
        for k = 1:length(pp.majorIter(i).minorIter)
           l = l + 1;
           for j = 1:N
                eci2lvlh    = lvlh2eci(pp.x_dProp(1:3,j), ...
                                       pp.x_dProp(4:6,j))';
                pMan(:,j) = eci2lvlh*pp.majorIter(i).minorIter(k).p(:,j);
            end
            plot3(pMan(1,:),pMan(2,:),pMan(3,:),'color',gradColors(l,ncolors))
            hold on
            leg = char([leg; string(['j = ',num2str(i) ', k = ',num2str(k)])]);
        end
    end
    legend(leg);
    xlabel('[m]')
    ylabel('[m]')
    zlabel('[m]')
    hold off
    grid on
end

%% Relative trajectory with ellipsoids
if p.relTrajEll && isfield(pp.majorIter(end).minorIter(end),"ellipsoid")
     relTrajEllipsoids(pp)
end

%% Relative trajectory with DeltaV
if p.relTrajDv
     relTrajDv(pp)
end

%% Acceleration
if p.acc
    figure()
    plot(t,pp.uRtn','LineWidth',2)
    hold on
    plot(t,pp.uNorm,'k-.','LineWidth',2)
    xlabel('Number of orbits [-]')
    ylabel('$u$ [-]')
    legend('$u_R$','$u_T$','$u_N$','$||u||$','Interpreter','latex')
%     title('Control action in RTN frame')
    grid on
    hold off
end

%% DeltaV
if p.dv
    figure()
    dv = pp.dv*1000;
    dv(abs(dv)< 1e-3) = 0;
    dvNorm = normOfVec(dv);
    if pp.uMax*pp.scaling(7) >= 1e-6                                 
        dv(dv==0) = nan;
        stem(t,dv(1,:),'LineWidth',2)
        hold on
        stem(t,dv(2,:),'LineWidth',2)
        stem(t,dv(3,:),'LineWidth',2)
        stem(t(dvNorm~=0),dvNorm(dvNorm~=0),'color','k','LineWidth',2)
        plot(t(dvNorm==0),dvNorm(dvNorm==0),'color','k')
    else
        plot(t,dv','LineWidth',2)
        hold on
        plot(t,dvNorm,'LineWidth',2,'color','k')
    end
    xlabel('Number of orbits [-]')
    ylabel('$\Delta v$ [mm/s]')
    legend('R','T','N','$|\cdot|$','interpreter','latex')
    grid on
    axis tight
    hold off
%     a = lvlh2rtn*pp.dv;
%     A = pp.DvNorm;
%     theta = atan2d((a(1,:).^2+a(2,:).^2),a(3,:));
%     phi = atan2d(a(1,:),a(2,:));
%     theta(theta==0) = nan; 
%     phi(phi==0) = nan; 
%     figure()
%     subplot(3,1,1)
%     if pp.uMax >= 1e-2;     A(A==0) = nan; stem(t,A)
%     else; plot(t,A)
%     end
%     xlabel('t/T [-]')
%     ylabel('\Delta V [m/s^2]')
%     grid on
%     axis tight
%     subplot(3,1,2)
%     if pp.uMax >= 1e-2; theta(theta==0) = nan; end; plot(theta,'*')
%     xlabel('t/T [-]')
%     ylabel('\theta [deg]')
%     grid on
%     axis tight
%     subplot(3,1,3)
%     if pp.uMax >= 1e-2; phi(phi==0) = nan; end; plot(phi,'*')
%     xlabel('t/T [-]')
%     ylabel('\phi [deg]')
%     grid on
%     axis tight
% end

if p.linCaConstr
    for i = pp.NCA0:pp.NCAf
        c(i) = dot(grad(:,i),pNew(:,i)-z(:,i)*pp.scaling(1));
    end
figure()
plot(pp.t/pp.T,c)
grid on
xlabel('Number of orbits []')
ylabel('$\nabla d_m^2\cdot(r-z)$')
end
end