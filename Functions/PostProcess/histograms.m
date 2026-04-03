function [] = histograms(var,edges)
% Build histogram variable
num = size(var,1);
for k = 1:num
    counts(k,:) = histcounts(var(k,:),edges);
end
counts = counts'/size(var,2)*100;

% Plot histograms
bar(edges(2:end),counts,'group','FaceColor','flat','EdgeColor','flat');
ylabel('\% of cases [-]')
grid on
box on
end