function plotSmd(t,pp)
%PlotSmd Plots the orginal and optimised SMD evolutions
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
M     = length(pp.secondary);
for j = 1:M
    figure()
    if ~strcmp(pp.obj,'miss_distance')
        semilogy(t,pp.ballisticSmd(:,j))
        hold on
        semilogy(t,pp.sqrMahaMan(:,j))
        semilogy(t,pp.sqrMahaLim(:,j),'k--')
        if pp.enableSmdGradConstraint
            semilogy(t,pp.sqrMahaLim(:,j)*1.1,'g--')
            legend({'Ballistic','Optimized','Limit',...
                    'Trigger \nabla(d_m^2) constraint'},'NumColumns',2)
    
        else
            legend('Ballistic','Optimized','Limit')
        end
        ylabel('$d_m^2$ [-]')
        axis tight
        ylim([min(pp.ballisticSmd(:,j))/10,max(pp.sqrMahaMan(:,j))*10])
    else
        plot(t,pp.ballisticSmd(:,j).^.5*pp.scaling(1))
        hold on
        plot(t,pp.sqrMahaMan(:,j).^.5*pp.scaling(1))
        plot(t,pp.sqrMahaLim(:,j).^.5*pp.scaling(1),'k--')
        ylabel('$d_{sep}$ [km]')
        legend('Ballistic','Optimized','$d_{sep}$ limit','interpreter','latex')
        axis tight
        ylim([min(pp.ballisticSmd(:,j).^.5*pp.scaling(1))/10, ...
              max(pp.sqrMahaMan(:,j).^.5*pp.scaling(1))*1.2])
    end
    xlabel('Number of orbits [-]')
    grid on
    hold off
end
end