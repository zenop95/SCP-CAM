function plotLimsEvolution(pp)

    color1 = [0 0.4470 0.7410]  ;
    color2 = [0.4660 0.6740 0.1880];
    color3 = [0.8500 0.3250 0.0980];
    colors = color1 + 1/4*(color2-color1).*repmat([0:4,4*ones(1,5)],3,1)' + ...
                      1/5*(color3-color2).*repmat([zeros(1,4),0:5],3,1)';
    lim = nan(length(pp.majorIter),length(pp.secondary));
    c = [];
    s = [];
    for kk = 1:pp.n_conj
        for jj = 1:pp.gmmOrder
            k = (kk-1)*pp.gmmOrder + jj; 
            for j = 1:length(pp.majorIter)
                lim(j,k) = pp.majorIter(j).lims(k);
            end
%             if pp.gmmOrder > 1; 
            c = [c;num2str(jj)]; s = [s;num2str(kk)]; 
%             end
        end
    end
    c = string(c); s = string(s);
     if pp.gmmOrder > 1 && pp.fastEncounter
        lim(:,ceil(pp.gmmOrder/2)) = sum(lim(:,1:pp.gmmOrder),2);
        lim(:,[1:ceil(pp.gmmOrder/2)-1,ceil(pp.gmmOrder/2)+1:pp.gmmOrder]) = [];
        c(1:pp.gmmOrder) = []; c = ["";c];
        s(1:pp.gmmOrder) = []; s = ["1";s];
    else
        s = string(num2str([1:pp.n_conj]'));
    end
    figure()
    for j = 1:size(lim,2)
        if pp.gmmOrder > 1 && pp.fastEncounter
            p = [char(c(j)), char(s(j))]; 
        elseif pp.gmmOrder == 1 && pp.fastEncounter
            p = char(s(j));
        else; p = char(c(j));
        end
        semilogy(lim(:,j),'o-','color',colors(j,:),'LineWidth',2, ...
                    'DisplayName',['$\bar{P}_{C,',p,'}$'])
        hold on
    end

    semilogy([1,length(pp.majorIter)],pp.lim*[1,1],'k--','DisplayName','$\bar{P}_{TC}$')
    if length(pp.majorIter)<7; xticks(1:length(pp.majorIter)); end
    leg = legend;
    set(leg, 'Interpreter','latex')
    grid on
    xlabel('Major Iterations [-]')
    if pp.gmmOrder > 1; c = 'c'; else c = []; end
    if pp.fastEncounter
        ylabel(['$\bar{P}_{C,',c,'s}$'])
    else
        ylabel(['$\bar{P}_{IC,',c,'s}$ [-]'])
    end
    ylim([pp.lim*1e-3,pp.lim*1.1])
    hold off
    
end