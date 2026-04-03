function  ref = evalTraj(DAmaps,N,nVars,dx)

% extract the unperturbed trajectory 
ref = nan(nVars,N);

for i =1:N
    for j = 1:nVars
        % exact position and velocity of the primary at each time step
        ref(j,i) = eval_poly(DAmaps(j,i).C, DAmaps(j,i).E,dx);
    end
end