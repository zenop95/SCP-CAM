
function p = eval_poly(C,E,PoI)
% PoI vettore riga 

p = [];
for i = 1:size(PoI,1)
    
    if isempty(C)
        p = [p;0];
    else
        PoImx = ones(length(C),1)*PoI(i,:);
        p = [ p; sum(C.*prod((PoImx.^E),2)) ];
    end
    
end