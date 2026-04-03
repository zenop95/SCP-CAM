load('data')
B = table2array(data); clear data;

means = mean(B(:,9:14));
varis = var(B(:,9:14));
varis(2) = mean([varis(1),varis(3)]);
c11 = abs(normrnd(means(1),varis(1),1,1,10));
c22 = abs(normrnd(means(2),varis(2),1,1,10));
c33 = abs(normrnd(means(3),varis(3),1,1,10));
c12 = normrnd(means(4),varis(4),1,1,10);
c13 = normrnd(means(5),varis(5),1,1,10);
c23 = normrnd(means(6),varis(6),1,1,10);
cov = [c11 c12 c13; c12 c22 c23; c13 c23 c33];