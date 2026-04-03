function [A,b] = smdGradBound(nOpt,NCA0,NCAf,gradOn,m,gradLim,pp)
aNew = zeros(1,nOpt);
A = []; b = [];
for i = NCA0:NCAf
    if gradOn(i)
        aNew(pp.index.gradCone + m*(i-1)) = 1;                          
        aNew(pp.index.gradVc + m*(i-1))   = 1;                         
        A = [A; aNew];                                               
        aNew = zeros(1,nOpt);
        b = [b; gradLim(i)/pp.gradScale];
    end
end
A = sparse(A);
end