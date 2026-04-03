function pp = generateInitSerra(pp)

%% Covariances
c11 = 67.013620376984; c12 = 14.572096438565; c13 = 31.362985178843;  c14 = -0.00183732409246960;  c15 = 0.00397493999345520; c16 = 0.00775436576273830;
                       c22 = 3.2132921306927; c23 = 6.8233552264667;  c24 = -0.00039892190369618;  c25 = 0.00085775667714172; c26 = 0.00168828264084460;
                                              c33 = 14.728803481810;  c34 = -0.00085890275305298;  c35 = 0.00183498030328000; c36 = 0.00363890409140320;
                                                                      c44 = 5.262586931470100e-8;  c45 = -1.0901907556364e-7; c46 = -2.1622909141821e-7;  
                                                                                                   c55 = 2.38807999319590e-7; c56 = 4.53189781842900e-7;
                                                                                                                              c66 = 9.07514528582240e-7;   
p.C0 = [c11        c12        c13         c14      c15          c16;
        c12        c22        c23         c24      c25          c26;
        c13        c23        c33         c34      c35          c36;
        c14        c24        c34         c44      c45          c46;
        c15        c25        c35         c45      c55          c56;
        c16        c26        c36         c46      c56          c66]/1e6;

c11 = 67.014285755464; c12 = 14.572190820724; c13 = 31.363196175928;  c14 = -0.00183733343291710;  c15 = 0.00391496994321290; c16 = 0.00775442526685390;
                       c22 = 3.2133009783397; c23 = 6.8233774594117;  c24 = -0.00039892244761992;  c25 = 0.00085775960971626; c26 = 0.00168828971574030;
                                              c33 = 14.728855328372;  c34 = -0.00085890413424338;  c35 = 0.00183435791016330; c36 = 0.00363891961183050;
                                                                      c44 =  5.26259182176670e-8;  c45 = -1.0901938617321e-7; c46 = -2.1622971004667e-7;  
                                                                                                   c55 = 2.38809208123100e-7; c56 = 4.53186178994030e-7;
                                                                                                                              c66 = 9.07519296888050e-7;
C0s = [c11        c12        c13         c14      c15          c16;
       c12        c22        c23         c24      c25          c26;
       c13        c23        c33         c34      c35          c36;
       c14        c24        c34         c44      c45          c46;
       c15        c25        c35         c45      c55          c56;
       c16        c26        c36         c46      c56          c66]/1e6;

%% Absolute state of the primary object at TCA
mu       = 398600.4418;    % [km^3/s^2]
p.HBR    = 0.003;           % [km]
p.mass   = 260;            % [kg] mass
p.A_drag = 1;              % [m^2] drag surface area
p.Cd     = 2.2;            % [-] shape coefficient for drag
p.A_srp  = 1;              % [m^2] SRP surface area
p.Cr     = 1.31;           % [-] shape coefficient for SRP
p.cart0  = [-5532.7006575059 20132.6739581200  40010.548546273 ...
            -1.4509451284357 -0.3116085722286 -0.6713019125023]';
p.x0     = p.cart0;
coe      = cartesian2kepler(p.cart0);
T        = 2*pi/coe.n;         % [s] orbital period
p.a      = coe.a;

%% Absolute state of the secondary object at TCA
j = 1;
s(j).tca      = 1;         % [s] TCA of conjunction w.r.t. initial time t0 = 0
s(j).x0       = [-5532.6940174556 20132.6765073780  40010.553862016 ...
                 -1.4509465079251 -0.3116078572248 -0.6713005315149]';
s(j).relState   = [];                        % [km] [km/s] Relative cartesian state at TCA
s(j).C0         = C0s;             % [km^2] [km^2/s^2] Covariance at TCA
s(j).HBR        = p.HBR + 0.003;  % [km]
s(j).mass       = 200;          % [kg] mass
s(j).A_drag     = 1;            % [m^2] drag surface area
s(j).Cd         = 2.2;          % [-] shape coefficient for drag
s(j).A_srp      = 1;            % [m^2] SRP surface area
s(j).Cr         = 1.31;         % [-] shape coefficient for SRP
s(j).x          = [];         
s(j).covariance = [];         
s(j).w          = 1;
s(j).cdm        = true;
s(j).ang        = false;
    
pp.fastEncounter = false;                                          
pp.primary       = p;
pp.secondary     = s;
pp.T             = T;
pp.initCovRtn    = true;

end