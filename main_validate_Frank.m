% 0) Path init / toolboxes
initializePath();
n           = 2170;
DvTot       = nan(n,1);
poc         = nan(n,1);
md          = nan(n,1);
t_thrust    = nan(n,1);

%%
for i = 1:2170
    try
        fprintf('\n============================================================\n');
        fprintf('Running scenario %d / %d \n', i, n);
        fprintf('============================================================\n');
        
        t_thrust(i) = 0;
        pp          = simProperties(i);
        pp          = runSCP(pp);
        t           = pp.t;                   % [-] (1,1) Time step
        dt          = diff(t); 
        dt(end+1)   = dt(end);
        dv          = pp.majorIter(end).uRtn.*repmat(dt',3,1)*pp.uMax;                        % [m/s] (3,N) Node-wise Delta V components in LVLH
        DvNorm      = normOfVec(dv);
        DvTot(i)    = sum(DvNorm)*pp.scaling(4)*1000;
        poc(i)      = pp.majorIter(end).pc;
        vp          = pp.majorIter.absV;
        rp          = pp.majorIter.absP;
        e2b         = eci2Bplane(vp(:,end),pp.secondary.cart(4:6,end));
        e2b         = e2b([1 3],:);
        md(i)       = norm(e2b*(rp(:,end)-pp.secondary.cart(1:3,end)))*pp.scaling(1);
        for j = 1:size(dv,2)
            if abs(1 - pp.uMax*dt(j)/DvNorm(:,j)) < 0.1
                t_thrust(i) = t_thrust(i)+dt(j)*pp.Tsc;
            else
                break;
            end
        end
    catch ME
        warning('Scenario %d failed: %s', i, ME.message);
    end
end
%%
save('FO_comparison_MD.mat', 'poc', 'DvTot', 't_thrust');
