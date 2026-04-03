function linear = extractLinear(DA)
%extractLinear Extracts the Linear Part of the DA polynomials DA.
for i = 1:length(DA)    
    % order of the DAMaps row (0 for the first, 1 for 2:10,)
    order = sum(DA(i).E,2);
    % take columns of DA correspondent to order 1
    lin = DA(i).E(order==1,:);
    coe = DA(i).C(order==1);
    linear(i,:) = sum(coe.*lin,1);
end