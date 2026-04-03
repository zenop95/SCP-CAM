function showGmm(mean,cov,ind)
figure
for j = 1:size(mean,3)
    showEllipsoid(mean(:,ind,j),cov(1:3,1:3,ind,j))
hold on
end
hold off
end