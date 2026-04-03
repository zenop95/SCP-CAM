function savePlot(name, path, type, res)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ax = gca; 
ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.ZTickMode = 'manual';
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPositionMode = 'manual';

print(strcat(fullfile(path,name),'.png'),type,res);
savefig(strcat(fullfile(path,name),'.fig'));
end

