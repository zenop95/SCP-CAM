function [] = absTrajDv(pp)
%RELTRAJDV Summary of this function goes here
%   Detailed explanation goes here\
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

N          = pp.N;
pNew       = pp.majorIter(end).absP;
for i = 1:N
    r2e     = rtn2eci(pNew(:,i),pp.majorIter(end).absV(:,i));
    dv(:,i) = r2e*pp.u(:,i);
end
% dv         = dv/pp.uMax*5; % scale the dv vector for better visualization

% figure()
plot3(pNew(1,:),pNew(2,:),pNew(3,:))
hold on
xlabel('x_{ECI} [-]')
ylabel('y_{ECI} [-]')
zlabel('z_{ECI} [-]')
grid on
for i = 1:N
    if norm(dv(:,i))
        quiver3(pNew(1,i),pNew(2,i),pNew(3,i),dv(1,i),dv(2,i),dv(3,i),'r','AutoScale','on')
    end
end
plot3(pNew(1,1),pNew(2,1),pNew(3,1),'o')
plot3(pNew(1,pp.N_back + 1),pNew(2,pp.N_back + 1),pNew(3,pp.N_back + 1),'*')
hold off

end

