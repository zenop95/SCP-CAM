function [nu_max,rho,Ju,nonLinCost,refute] = updateXiMax(pp,majorIter,nu_max,linTraj,iter)
N       = pp.N;
rho1    = 0.1;
rho2    = 0.2;
rho3    = 0.4;
rho4    = 0.8;
beta1   = 1/4;
beta2   = 4;
% nonlinear keplerian propagation
% nonLinX(:,1) = linTraj(:,1);
% for i = 1:pp.N-1
%     uEci = rtn2eci(linTraj(1:3,i),linTraj(4:6,i))*majorIter.u(:,i)*pp.uMax;
%     nonLinX(:,i+1) = propKepOde(linTraj(:,i),uEci,pp.t(i+1)-pp.t(i),1);
% end
u       = reshape(majorIter.u,3*N,1);
fid = fopen('initial_state.dat', 'w');
fprintf(fid, '%2i\n',     N);
fprintf(fid, '%40.12f\n', pp.t);
fprintf(fid, '%40.12f\n', pp.et);
fprintf(fid, '%2i\n',     1);
fprintf(fid, '%40.12f\n', pp.uMax);
fprintf(fid, '%40.12f\n', pp.scaling(1));
fprintf(fid, '%2i\n',     pp.flagRtn); %Flag for type of propagation
fprintf(fid, '%2i\n',     2); %Flag for type of propagation
for i = 1:6*N
     fprintf(fid, '%40.12f\n', linTraj(i));       
end
for i = 1:3*N
    fprintf(fid, '%40.12f\n', u(i));
end
fclose(fid);
aidaInit(pp,'primary');
!wsl ./CppExec/propAida
nonLinX      = reshape(load("constCart.dat"),6,N);
pp.timeSubtr = pp.timeSubtr + load("timeOut.dat");

linX    = reshape(linTraj,6,N);
nonLinVc  = max(normOfVec(nonLinX-linX))*pp.vcWeight;
Jvc = max(majorIter.minorIter(end).virtual)*pp.vcWeight; 
Jt    = 0;
nonLinTarg = 0;
if pp.enableSkTarget
    Jt    = pp.targWeight*norm(majorIter.targSlack);
    nonLinTarg = pp.targWeight*norm(nonLinX(:,end) - pp.skTarget);
end

%% Update
Ju         = sum(normOfVec(majorIter.u))*pp.ctrlWeight;
nonLinCost = Ju + nonLinVc + nonLinTarg;
linCost    = Ju + Jvc + Jt;
if iter > 1
    DeltaJ  = pp.majorIter(iter-1).J - nonLinCost;
    DeltaL  = pp.majorIter(iter-1).J - linCost;
else
    refute = false;
    rho    = 1;
    return
end

% if DeltaJ>-1e-6 && DeltaJ<0; DeltaJ = -DeltaJ; end
rho     = DeltaJ/DeltaL;
refute  = false;
if nu_max == pp.nu_m
    refute = false;
    return
end

beta = computeBeta(rho,rho1,rho3,rho4,beta1,beta2);
nu_max = min(pp.nu_M,max(pp.nu_m,nu_max*beta));

if  rho <= rho2; refute = true; end

end

function beta = computeBeta(rho,rho1,rho3,rho4,beta1,beta2)
    
    alpha1 = (1-beta1)/(rho3-rho1);
    alpha2 = (beta2-1)/(1-rho4);

    if  rho <= rho1
        beta = beta1;
    elseif rho >= rho1 && rho < rho3
        beta = beta1 + alpha1*(rho-rho1);
    elseif rho >= rho4 && rho < 1
        beta = 1 + alpha2*(rho-rho4);
    elseif rho >= 1
        beta = beta2;
    else
        beta = 1;
    end
end