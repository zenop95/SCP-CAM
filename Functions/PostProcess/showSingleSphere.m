function [] = showSingleSphere(pp)
%relTrajEllipsoids plots the relative trajectories with ellipsoid
%visualization
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
lastIter   = pp.majorIter(end);
M          = length(pp.secondary);
color1 = [0 0.4470 0.7410]  ;
color2 = [0.4660 0.6740 0.1880];
colors = color1 + 1/M*(color2-color1).*repmat(0:M,3,1)';
    
figure()
hold on
for j = 1:length(pp.secondary)
    pNew       = pp.validationTraj(1:3,:,j);
    pOld       = pp.ballisticRelTraj(1:3,:,j);
    rN         = nan(3,pp.N);
    rO         = nan(3,pp.N);
    sqLims     = pp.majorIter(end).sqrMahaLim(:,j);
    P          = lastIter.P(:,:,:,j);
    if strcmpi(pp.obj,'miss_distance')
        ind = normOfVec(pp.ballisticRelTraj(1:3,:)) < pp.lim;
    else
        ind = pp.ballisticIpc(:,j) > pp.lim;
    end
    for i = 62:90
        if ind(i)
            smdLim     = sqLims(i);
            PEci       = P(:,:,i);
            [semiax,cov2e] = defineEllipsoid(PEci,smdLim);
            rN(:,i)  = cov2e'*pNew(:,i)./semiax;
            rO(:,i)  = cov2e'*pOld(:,i)./semiax;
        end
    end
    for i = 64:90
        ang = 0.001;
        ax = [1 0 0]';
        kk = 0;
        while norm(rN(:,i)-rN(:,i-1)) > norm(rN(:,i-1)-rN(:,i-2))*1.1
            rN(:,i) = axang2dcm(ax,ang)*rN(:,i);
            kk = kk + 1;
            if kk > pi/0.001
                ax = [0 1 0]';
            end
        end
%         while norm(rO(:,i)-rO(:,i-1)) > norm(rO(:,i-1)-rO(:,i-2))*1.1
%             rO(:,i) = axang2dcm(ax,ang)*rO(:,i);
%             kk = kk + 1;
%             if kk > pi/0.001
%                 ax = [0 1 0]';
%             end
%         end
    end
    plot3(rO(1,ind),rO(2,ind),rO(3,ind),'-','marker','diamond','Color', ...
        colors(j,:),'LineWidth',1,'DisplayName',['Ball c=',num2str(j)])
    plot3(rN(1,ind),rN(2,ind),rN(3,ind),'-o','Color', ...
        colors(j,:),'LineWidth',1.5,'DisplayName',['Opt c=',num2str(j)])
end
[ellX,ellY,ellZ] = ellipsoid(0,0,0,1,1,1,20);
surf(ellX,ellY,ellZ,'FaceAlpha',0.5,'LineStyle','none', ...
            'FaceColor',[0.7 0.7 0.7],'HandleVisibility','off');
plot3(0,0,0,'k*','LineWidth',2,'HandleVisibility','off')
xlabel('$\tilde{\xi}  $ [-]')
ylabel('$\tilde{\eta} $ [-]')
zlabel('$\tilde{\zeta}$ [-]')
grid on
view(3)
legend show
hold off
end