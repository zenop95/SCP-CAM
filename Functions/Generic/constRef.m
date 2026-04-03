function  ref = constRef(DAmaps,N,nVars)

% extract the unperturbed trajectory 
dx0NoMan  = zeros(1,size(DAmaps(1,1).E,2));
ref = nan(nVars,N);

for i =1:N
    for j = 1:nVars
        % exact position and velocity of the primary at each time step
        ref(j,i) = eval_poly(DAmaps(j,i).C, ...
            DAmaps(j,i).E,dx0NoMan);
    end
end