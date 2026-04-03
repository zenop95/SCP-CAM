function minorIter = checkCones(xNew,minorIter,pp)
% CheckCones checks that the optimization variables that define the cones 
% radiuses are close enough to the norm of the variables that are subjected
% to the cone constraint.
%
% INPUT: xNew      = [m*N] Solution of the optimization process
%        minorIter = [ ] Minor iteration structure.
%        pp        = [ ] Postprocess structure.
% 
% OUTPUT: miorIter = [ ] Minor iteration structure.
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
ind = pp.index;
% if any((xNew(ind.ctrlCone,:)-normOfVec(xNew(ind.ctrl,:)))>1e-5)
%     warning('The acceleration cone weight is too low, consider increasing it')
% end
minorIter.accCost = sum(xNew(ind.ctrlCone,:)*pp.ctrlWeight);
minorIter.vc   = xNew(ind.stateVc,:);
minorIter.vbCa = xNew(ind.caVc,:);
minorIter.vbGr = xNew(ind.gradVc,:);
minorIter.vbSk = xNew(ind.skVc,:);            
%     minorIter.virtual = xNew(ind.vcCone,:);
minorIter.virtual = normOfVec(xNew(ind.vc,:));
minorIter.vcCost = pp.vcWeight*sum(minorIter.virtual);
minorIter.trCost  = 0;
minorIter.gradCost  = 0;
minorIter.cost = minorIter.accCost + minorIter.trCost + ...
                 minorIter.vcCost + minorIter.gradCost;