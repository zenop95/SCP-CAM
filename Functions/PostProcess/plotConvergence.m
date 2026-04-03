function plotConvergence(pp)
%PlotConvergence plots the convergence evolution of major and minor
%iterations
% the SMD
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
p        = pp.figs;
lastIter = pp.majorIter(end);
if p.minorConvergence
    nMinorIter = length(lastIter.minorIter);
    dvTotMin   = nan(nMinorIter,1);
    errMin     = nan(nMinorIter,1);
    for j = 1:nMinorIter
        dvTotMin(j) = lastIter.minorIter(j).DvTot;
        errMin(j)   = lastIter.minorIter(j).err;
    end
    
    figure()
    plot(1:nMinorIter,dvTotMin)
    title('Minor Iterations \DeltaV profile')
    xlabel('Minor Iteration')
    ylabel('Total $\DeltaV$ [m/s]')
    grid on
    
    figure()
    semilogy(1:nMinorIter,errMin)
    title('Minor Iterations error profile')
    xlabel('Minor Iteration')
    ylabel('Iteration Error [m]')
    grid on
end

if p.majorConvergence
    nMajorIter = length(pp.majorIter);
    dvTotMaj   = nan(nMajorIter,1);
    errMaj     = nan(nMajorIter,1);
    for j = 1:nMajorIter
        dvTotMaj(j) = pp.majorIter(j).DvTot;
        errMaj(j)   = pp.majorIter(j).err;
        pc(j)       = pp.majorIter(j).PcTot;
        nu(j)       = pp.majorIter(j).nu_max;
    end

    figure()
    plot(1:nMajorIter,dvTotMaj)
    xlabel('Major Iteration')
    ylabel('Total $\Delta V$ [m/s]')
    grid on
    
    figure()
    semilogy(1:nMajorIter,errMaj)
    xlabel('Major Iteration')
    ylabel('Iteration Error [m]')
    grid on

    figure()
    semilogy(1:nMajorIter,pc)
    xlabel('Major Iteration')
    ylabel('Total $P_{TC}$ [-]')
    grid on

    figure()
    semilogy(1:nMajorIter,nu)
    xlabel('Major Iteration')
    ylabel('$\nu_{max}$ [-]')
    grid on
end
end