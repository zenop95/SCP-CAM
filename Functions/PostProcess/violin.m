function A = violin(var,varargin)
% Violin Plot
if nargin == 1
    A = iosr.statistics.boxPlot(var);
else
    A = iosr.statistics.boxPlot(varargin{1},var);
end
A.mediancolor = 'b';
A.boxWidth = 0;
A.outliersize = nan;
A.showMean = true;
A.percentile = [50 50]; 
A.lineStyle = 'none';
A.meancolor = 'b';
A.showViolin = true;
A.showScatter = false;
A.violinAlpha = 0.1;
A.LineColor   = [0 0.4470 0.7410];
A.violinColor = [0 0.4470 0.7410];
grid on
box on
end