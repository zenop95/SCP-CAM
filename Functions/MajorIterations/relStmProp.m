function [relTraj,uP,pp] = ...
                      relStmProp(iter,state,cart,pp)
% relStmProp propagates variables for major iteration
%
% INPUT: 
% 
% OUTPUT: 
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
N = pp.N;
if iter == 1 
    pp.originalState = state;                                
end
% Relative trajectory
relTraj = cart - pp.secondary.cart;

%% Check controllability node by node
% for i = 2:N
%     Co = ctrb(STM(:,:,i),dynMaps(:,1:3,i));
%     unco = 6 - rank(Co);
%     if unco > 2
%         warning(['the system is uncontrollable in node ', num2str(i)]);
%     end
% end

%% Homotopy
if pp.enableHomotopy
    if iter == 1
        uP = zeros(1,N);
    else
        uP = pp.majorIter(iter-1).uP;
    end
else
    uP = [];
end
end