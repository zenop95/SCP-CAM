function showEllipsoids(pp)
N = pp.N;
lastIter = pp.majorIter(1);
if isfield("InitEllipsoids",lastIter)
    for i = 1:N
        figure(1)
        plotEllipsoidInit(lastIter.InitEllipsoids(i))
        title('Initial expanded point')
    end
end
for i = 2:N
    plotSingleEllipsoid(lastIter,i);
end
plotSingleEllipsoid(lastIter,N);
for i = 1:N
    for j = 1:length(lastIter.minorIter)
        if isfield("ellipsoids",lastIter.minorIter(j))
            figure(j)
            plotEllipsoidPoints(lastIter.minorIter(j).ellipsoids(i), ...
                                    lastIter.minorIter(j).ellipsoids(i+1))
            title(['Minor iteration ',num2str(j)])
            hold on
        end
    end
end
hold off
end