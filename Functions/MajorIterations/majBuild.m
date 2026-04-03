function majorIter = majBuild(majorIter,t,nu_max,prob, ...
    dynMaps,llaMaps,cartMaps,probMaps,state,lla,constPos,C,P,smdLim, ...
    highRisk,gradLim,gradOn,xi_dyn,xi_ca,pp)
% MajBuild builds the structure of the current major iteration.
%
% INPUT: majorIter = [struct]    Structure containing output major iteration
%        ll_const  = [-] (6,1)   Constant part of the geodetic coordinates
%        P         = [-] (6,6,N) Positional covariance matrices
%        smdLim    = [-] (N,1)   Node-wise SMD limit   
%        STM       = [-] (6,6,N) STM of each node
%        nu        = [-] (6,N)   Nonlinearity parameter nu
%        pp        = [struct]    Postprocess structure
%        
% OUTPUT: majorIter = [struct]  Structure containing output major iteration
%         
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
majorIter.err        = 1;
majorIter.J          = 1;
majorIter.linCost    = 1;
majorIter.x          = majorIter.minorIter(end).x;
majorIter.absP       = majorIter.minorIter(end).p;
majorIter.absV       = majorIter.minorIter(end).v;
M                    = length(pp.secondary);
majorIter.relP = majorIter.minorIter(end).r_rel;
for j = 1:M
    majorIter.relV(:,:,j) = majorIter.absV - pp.secondary(j).cart(4:6,:);
end
majorIter.u          = majorIter.minorIter(end).u;
majorIter.uNorm      = majorIter.minorIter(end).uNorm;
for i = 1:pp.N
    r2e = rtn2eci(majorIter.absP(:,i),majorIter.absV(:,i));
    if pp.flagRtn
        majorIter.uRtn(:,i) = majorIter.u(:,i);
        majorIter.uEci(:,i) = r2e*majorIter.u(:,i);
    else
        majorIter.uEci(:,i) = majorIter.u(:,i);
        majorIter.uRtn(:,i) = r2e'*majorIter.u(:,i);
    end
end
majorIter.uP         = majorIter.minorIter(end).uP;
majorIter.xi         = majorIter.minorIter(end).xi;
majorIter.DvNorm     = majorIter.uNorm*pp.dt*pp.uMax*pp.scaling(4);
majorIter.DvTot      = sum(majorIter.DvNorm);
majorIter.sqrMahaLim = smdLim;

if pp.fastEncounter
    [PoC,smd]  = computePoC(majorIter.relP,P,pp);
    ipc        = nan;
else
    [ipc,smd]  = computeIpc(majorIter.relP,P,pp);
    PoC        = nan;
end
majorIter.sqrMaha    = smd;
majorIter.ipc        = ipc;
majorIter.ipcTot     = ipcTot(ipc);
majorIter.pc         = PoC;
majorIter.PcTot      = PoCTot(PoC);
majorIter.covPrim    = C;
majorIter.P          = P;
majorIter.deviation  = normOfVec(majorIter.xi./nu_max);
majorIter.nu_max     = nu_max;
majorIter.gradLim    = gradLim;
majorIter.lims       = pp.lims;
majorIter.t          = t;
majorIter.xi_dyn     = xi_dyn;
majorIter.xi_ca      = xi_ca;
majorIter.gradLim    = gradLim;
majorIter.dynMaps    = dynMaps;
majorIter.llaMaps    = llaMaps;
majorIter.cartMaps   = cartMaps;
majorIter.probMaps   = probMaps;
majorIter.state      = state;
majorIter.uP         = majorIter.minorIter(end).uP;
majorIter.lla        = lla;
majorIter.constPos   = constPos;
majorIter.highRisk   = highRisk;
majorIter.gradOn     = gradOn;
majorIter.prob       = prob;
if pp.enableSkTarget
    majorIter.targSlack  = majorIter.minorIter(end).slack;
end
end