function plotLlaEcef(pp,t)
%PlotLlaEcef Plots the Latitude, longitude, altitude and ECEF variables
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
p         = pp.figs;
firstIter = pp.majorIter(1);
N         = pp.N;
llaOld = firstIter.lla; llaOld(1:2,:) = rad2deg(llaOld(1:2,:));
llaNew = pp.lla;        llaNew(1:2,:) = rad2deg(llaNew(1:2,:));
llaSK = repmat([0; rad2deg(pp.nomLon); llaOld(3,1)],[1,N]);

%% plots
if p.ecef 
    [ecefSKup(1),ecefSKup(2),ecefSKup(3)] = lla2ecef(deg2rad(llaSK(1,1) + ...
                pp.skDev(1)),deg2rad(llaSK(2,1)+pp.skDev(2)),llaSK(3,1)+pp.skDev(3));
    [ecefSKlo(1),ecefSKlo(2),ecefSKlo(3)] = lla2ecef(deg2rad(llaSK(1,1) - ...
                pp.skDev(1)),deg2rad(llaSK(2,1)-pp.skDev(2)),llaSK(3,1)-pp.skDev(3));
    for i = 1:N
        [ecefOld(1,i),ecefOld(2,i),ecefOld(3,i)] = ...
           lla2ecef(deg2rad(llaOld(1,i)),deg2rad(llaOld(2,i)),llaOld(3,i));
        [ecefNew(1,i),ecefNew(2,i),ecefNew(3,i)] = ...
           lla2ecef(deg2rad(llaNew(1,i)),deg2rad(llaNew(2,i)),llaNew(3,i));
    end
end
if p.lla
    if strcmpi(pp.orbit,'geo') &&  pp.enableSkTarget
        t_aft     = pp.t_aft'/pp.T;
        t_augm   = [t;t_aft];
        figure()
        subplot(2,1,1)
        hold on
        plot(t_augm,[llaOld(1,:), pp.lat_aftNoMan])
        plot(t_augm,[llaNew(1,:), pp.lat_aftMan])
        plot(t_augm,llaSK(1,1) + pp.skDev(1)*ones(length(t_augm),1),'k--')
        plot(t_augm,llaSK(1,1) - pp.skDev(1)*ones(length(t_augm),1),'k--')
        xlabel('Number of orbits [-]')
        ylabel('Latitude [deg]')
        axis tight
        legend('Original','Optimised','Limit')
        grid on
        hold off
        subplot(2,1,2)
        hold on
        plot(t_augm,[llaOld(2,:), pp.lon_aftNoMan])
        plot(t_augm,[llaNew(2,:), pp.lon_aftMan])
        plot(t_augm,llaSK(2,1) + pp.skDev(2)*ones(length(t_augm),1),'k--')
        plot(t_augm,llaSK(2,1) -pp.skDev(2)*ones(length(t_augm),1),'k--')
        xlabel('Number of orbits [-]')
        ylabel('Longitude [deg]')
        axis tight
        legend('Original','Optimised','Limit')
        grid on
        hold off
    else
       figure()
        subplot(2,1,1)
        hold on
        plot(t,llaOld(1,:))
        plot(t,llaNew(1,:))
        plot(t,llaSK(1,1) + pp.skDev(1)*ones(length(t),1),'k--')
        plot(t,llaSK(1,1) - pp.skDev(1)*ones(length(t),1),'k--')
        xlabel('Number of orbits [-]')
        ylabel('Latitude [deg]')
        axis tight
        legend('Original','Optimised','Limit')
        grid on
        hold off
        subplot(2,1,2)
        hold on
        plot(t,llaOld(2,:))
        plot(t,llaNew(2,:))
        plot(t,llaSK(2,1) + pp.skDev(2)*ones(length(t),1),'k--')
        plot(t,llaSK(2,1) -pp.skDev(2)*ones(length(t),1),'k--')
        xlabel('Number of orbits [-]')
        ylabel('Longitude [deg]')
        axis tight
    %     title('Deviations from initial longitude')
        legend('Original','Optimised','Limit')
        grid on
        hold off
    end
    figure()
subplot(1,2,1)
    plot([llaOld(2,:), pp.lon_aftNoMan],[llaOld(1,:), pp.lat_aftNoMan])
    hold on
    plot(linspace(llaSK(2,1)-pp.skDev(2),llaSK(2,1)+pp.skDev(2),N),(llaSK(1,1)+pp.skDev(1))*ones(N,1),'k--')
    plot(linspace(llaSK(2,1)-pp.skDev(2),llaSK(2,1)+pp.skDev(2),N),(llaSK(1,1)-pp.skDev(1))*ones(N,1),'k--')
    plot((llaSK(2,1)+pp.skDev(2))*ones(N,1),linspace(llaSK(1,1)-pp.skDev(1),llaSK(1,1)+pp.skDev(1),N),'k--')
    plot((llaSK(2,1)-pp.skDev(2))*ones(N,1),linspace(llaSK(1,1)-pp.skDev(1),llaSK(1,1)+pp.skDev(1),N),'k--')
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')
    grid on
    hold off
subplot(1,2,2)
    plot([llaNew(2,:), pp.lon_aftMan],[llaNew(1,:), pp.lat_aftMan])
    hold on
    plot(linspace(llaSK(2,1)-pp.skDev(2),llaSK(2,1)+pp.skDev(2),N),(llaSK(1,1)+pp.skDev(1))*ones(N,1),'k--')
    plot(linspace(llaSK(2,1)-pp.skDev(2),llaSK(2,1)+pp.skDev(2),N),(llaSK(1,1)-pp.skDev(1))*ones(N,1),'k--')
    plot((llaSK(2,1)+pp.skDev(2))*ones(N,1),linspace(llaSK(1,1)-pp.skDev(1),llaSK(1,1)+pp.skDev(1),N),'k--')
    plot((llaSK(2,1)-pp.skDev(2))*ones(N,1),linspace(llaSK(1,1)-pp.skDev(1),llaSK(1,1)+pp.skDev(1),N),'k--')
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')
    grid on
    hold off
end

if p.ecef
    figure()
    plot3(ecefOld(1,:),ecefOld(2,:),ecefOld(3,:))
    hold on
    plot3(ecefNew(1,:),ecefNew(2,:),ecefNew(3,:))
    plot3dbox(ecefSKlo,ecefSKup)
    legend('Original','Optimised','Keep-in box')
    xlabel('[m]')
    ylabel('[m]')
    zlabel('[m]')
    title('ECEF trajectory')
    hold off
    grid on
end
end