function pp = findFiringNodes(pp)

pp.canFire = true(pp.N,1);
if ~pp.canAlwaysFire
    window = pp.fireWindow;
    start  = utc2et(window(1,:))/pp.Tsc;
    ending = utc2et(window(2,:))/pp.Tsc;
    toTca1  = start  - pp.et;
    toTca2  = ending - pp.et;
    Nstart  = round(toTca1/pp.dt) + 1;
    Nending = round(toTca2/pp.dt) + 1;
    pp.canFire(1:Nstart-1) = false;
    pp.canFire(Nending+1:end) = false;
end
end
