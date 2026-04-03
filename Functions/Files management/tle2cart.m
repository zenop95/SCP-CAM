function cart = tle2cart(tle)
    mu = 398600;
    s = size(tle,1);
    inc = nan(s,1);
    Om = nan(s,1);
    ecc = nan(s,1);
    om = nan(s,1);
    meanAn = nan(s,1);
    trueAn = nan(s,1);
    n = nan(s,1);
    a = nan(s,1);
    cart = nan(s,6);
    for i = 1:s
        ch         = tle(i);
        ch         = ch{:};
        inc(i)     = deg2rad(str2num(ch(1,80:86)));
        Om(i)      = deg2rad(str2num(ch(1,88:95)));
        ecc(i)     = str2num(['.',ch(1,97:103)]);
        om(i)      = deg2rad(str2num(ch(1,105:112)));
        meanAn(i)  = deg2rad(str2num(ch(1,114:121)));
        trueAn(i)  = mean2trueAnomaly(meanAn(i),ecc(i));
        n(i)       = str2num(ch(1,123:133))*2*pi/86400;
        a(i)       = (mu/n(i)^2)^(1/3);
        cart(i,:)  = kepler2cartesian(a(i),ecc(i),Om(i),inc(i),om(i),trueAn(i),mu);
    end
end