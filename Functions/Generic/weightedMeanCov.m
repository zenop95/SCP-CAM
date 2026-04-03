function [mu,cov] = weightedMeanCov(xx,w)
w = w/sum(w);
mu = sum(w.*xx);
coeff = 1/(1-w'*w);

for indi = 1:size(mu,2)
    for indj = 1:size(mu,2)
        row = w.*(xx(:,indi)-mu(:,indi));
        col = xx(:,indj)-mu(:,indj);
        cov(indi,indj) = coeff*row'*col;
    end
end
end

