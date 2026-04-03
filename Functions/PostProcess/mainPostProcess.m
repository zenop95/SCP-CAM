function [] = mainPostProcess(pp)
%MainPostProcess calls the required post process functions
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
if any(pp.majorIter(end).minorIter(end).virtual>1e-8)
    warning('The original non-convex problem is infeasible')
end
p = pp.figs;
t = pp.t/pp.T;

%% Display results
dispResults(pp);

%% Ellipsoids
if p.ellipsoids; showEllipsoids(pp); end
if p.singleSphere; showSingleSphere(pp); end

%% Convergence
plotConvergence(pp);

%% Trajectories
plotTrajectories(t,pp);

%% B plane for short-term
if pp.fastEncounter 
if pp.n_conj > 1
    showSingleCircle(pp); 
else
    showEllipseBplane(pp);
end
end

%% SMD gradient
if p.grad && pp.enableSmdGradConstraint && ~strcmp(pp.obj,'miss_distance'); plotSmdGrad(t,pp); end

%% ECI ECEF and LLA of initial and final trajectory
if p.lla || p.ecef; plotLlaEcef(pp,t); end
if p.alt
    plot(t,(normOfVec(pp.ballisticTraj(1:3,:))-1)*pp.scaling(1))
    hold on
    grid on
    plot(t,(normOfVec(pp.validationAbsTraj(1:3,:))-1)*pp.scaling(1))
    xlabel('Number of orbits [-]')
    ylabel('Altitude [km]')
    legend('Ballistic','Optimized')
    plot([t(1),t(end)],pp.altLim(1)*pp.scaling(1)*[1,1],'k--')
    plot([t(1),t(end)],pp.altLim(2)*pp.scaling(2)*[1,1],'k--')
    hold off
end
%% Collision indicators
% if p.sqMaha && p.ipc
%     figure(); subplot(2,1,1); plotSmd(t,pp); subplot(2,1,2); plotIpc(t,pp); 
% else
if p.sqMaha; plotSmd(t,pp); end
if p.ipc 
    if strcmp(pp.obj,'miss_distance') 
        plotMd(t,pp);
    else
        plotIpc(t,pp);
    end
% end

if p.nli
    figure()
    grid on
    hold off
    subplot(3,1,1)
    for i = 1:length(pp.majorIter)
        semilogy(t,normOfVec(pp.majorIter(i).xi))
        hold on
    end
    xlabel('$t/T$ [-]')
    ylabel('$\xi$ [-]')
    grid on
    hold off

    subplot(3,1,2)
    for i = 1:length(pp.majorIter)
        semilogy(t,normOfVec(pp.majorIter(i).nu))
        hold on
    end
    xlabel('$t/T$ [-]')
    ylabel('$\nu$ [-]')
    grid on
    hold off

    subplot(3,1,3)
    for i = 1:length(pp.majorIter)
        leg = ['iter ', num2str(i)];
        semilogy(t,normOfVec(pp.majorIter(i).xi./pp.majorIter(i).nu),...
                                                'DisplayName',leg)
        hold on
    end
    legend show
    xlabel('$t/T$ [-]')
    ylabel('$|\delta x|$ [-]')
    grid on
    hold off
end
if p.PoCTime
    figure()
    pocSingle = pp.PoCChan;
    pocSingleBall = pp.ballisticPoc;
    pocSummed = [1e-30; nan(pp.n_conj,1)];
    pocBall  = [1e-30; nan(pp.n_conj,1)];
    for j = 2:pp.n_conj+1
        pocSummed(j) = PoCTot(pocSingle(1:j-1));
        pocBall(j) = PoCTot(pocSingleBall(1:j-1));
    end
    semilogy(pp.t(pp.NCA)/pp.T,pp.ballisticPoc,'diamond','LineWidth',2,'color','k')
    hold on
    in = interp1([pp.t(1); pp.t(pp.NCA); pp.t(end)],[pocBall; pocBall(end)],pp.t,'previous');
    semilogy(pp.t/pp.T,in,'k','LineWidth',2)
    semilogy(pp.t(pp.NCA)/pp.T,pocSingle,'o','LineWidth',2,'color',[0 0.4470 0.7410])
    in = interp1([pp.t(1); pp.t(pp.NCA); pp.t(end)],[pocSummed; pocSummed(end)],pp.t,'previous');
    semilogy(pp.t/pp.T,in,'LineWidth',2,'color',[0 0.4470 0.7410])
    semilogy([pp.t(1)/pp.T,pp.t(end)/pp.T],pp.lim*ones(1,2),'k--');
    grid on
    legend('$P_C$ before the maneuver','$P_C$ after the maneuver', '$P_{TC}$','Interpreter','latex')
    ylim([pp.lim/1e4,1])
    xlabel('$t/T$ [-]')
    ylabel('$P_{TC}$ [-]')
    hold off
end

%% Angle from ellipse semi-major axis
if p.gamma
    ell = pp.majorIter(end).minorIter(end).ellipsoid;
    r   = pp.validationTraj(1:3,:);
    aEci = zeros(3,length(pp.NCA0:pp.NCAf));
    bEci = zeros(3,length(pp.NCA0:pp.NCAf));
    cEci = zeros(3,length(pp.NCA0:pp.NCAf));
    gamma1 = nan(length(pp.NCA0:pp.NCAf),1);
    gamma2 = nan(length(pp.NCA0:pp.NCAf),1);
    for i = pp.NCA0:pp.NCAf
        if pp.majorIter(end).sqrMahaLim(i)>1e-4
            aEci(:,i) = ell(i).cov2CW*[1;0;0];
            cEci(:,i) = ell(i).cov2CW*[0;0;1];
            bEci(:,i) = ell(i).cov2CW*[0;1;0];
            gamma1(i) = acosd(dot(r(:,i),aEci(:,i))/norm(r));
            gamma2(i) = acosd(dot(r(:,i),bEci(:,i))/norm(r));
        end
    end
    figure()
    plot(pp.t(2:end)/pp.T,gamma1(2:end),'*')
    hold on
    plot(pp.t(2:end)/pp.T,gamma2(2:end),'*')
    legend('angle from major','angle from minor');
    xlabel('t/T []')
    ylabel('[deg]')
    grid on
end

if p.limsEvolution && (pp.gmmOrder > 1 || pp.n_conj > 1)
    plotLimsEvolution(pp);
end

end