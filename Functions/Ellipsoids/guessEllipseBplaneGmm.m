function pointEci = guessEllipseBplaneGmm(rEci,smdLim,P,conj_n,pp)
% pOnEllipsoidLine draws a segment from a point in the interior of the 
% ellipsoid to the ellipsoid border and finds the point of intersection 
% between this line and the ellipsoid border.
%
% INPUT: rEci   = [m]   relative distance between the two bodies expressed in 
%                       the ECI frame centered on the center of the ellipsoid.                    
%        smdLim = [m^2] objective value of the squared Mahalanobis dist.
%        PEci   = [m^2] covariance in ECI
%        pp     = []    postprocess structure
% OUTPUT: pointEci  = [m] 3x1 vector of the position of the point on the 
%                         ellipsoid border expressed in CW coordinates.
%         ell       = [ ] Ellipsoid structure.
%
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
%pass in 2d to solve B-plane problem
ord        = pp.gmmOrder;
n          = 1e4;
centerGmm  = (conj_n-1)*ord + ceil(ord/2);
center     = pp.secondary(centerGmm).cart(1:3,pp.NCA(centerGmm));
a          = nan(2,ord);
offEci     = nan(3,ord);
ell        = nan(3,n,ord);
ellEci     = nan(3,n);
alpha      = nan(ord,1);
eta1       = pp.e2b(2,:,centerGmm)';
for j = 1:ord
    jj         = j+(conj_n-1)*ord;
    i          = pp.NCA(jj);
    PEci       = P(:,:,i,jj);
    e2b        = pp.e2b(:,:,jj);   % B space conversion for the jth mixand
    e2b2       = e2b([1,3],:);     % B plane conversion for the jth mixand
    PB2d       = e2b2*PEci*e2b2';  % Covariance in the B plane
    offset     = pp.secondary(jj).cart(1:3,pp.NCA(jj)) - center; % ECI Vector from center of central mixand to center of jth mixand
    eta2       = pp.e2b(2,:,jj)';  % ECI Perpendicular direction to the jth B-plane
    d          = dot(offset,eta1)/dot(eta1,eta2)*eta2; % ECI Vector from center of jth mixand to central B-plane perpendicular to jth B-plane
%      d=zeros(3,1);
    offEci(:,j)  = offset - d;     % ECI Vector from center of central mixand to intersection between the normal direction of the jth B-plane and the central B-plane

    [semiaxes,cov2b] = defineEllipsoid(PB2d,smdLim(jj));
    a(:,j)     = semiaxes;
    t          = linspace(0,2*pi,n);
    cov2e      = e2b2'*cov2b;
    alpha(j)   = acos(cov2b(1,1)); if cov2b(2,1) < 0; alpha(j) = -alpha(j); end
    ellCov     = [semiaxes(1)*cos(t); semiaxes(2)*sin(t)];
    for k = 1:n
        ellEci(:,k) = cov2e*ellCov(:,k) + offEci(:,j);
    end
    ell(:,:,j) = ellEci;
end

% Build outer envelope of the three ellipses
for j = 1:ord
    jj    = j + (conj_n - 1)*ord;
    e2b  = pp.e2b([1 3],:,jj);
    for kk = [1:j - 1, j + 1:ord]
        k  = kk + (conj_n - 1)*ord;
        %rotated and translated ellipse equation
        xc   = offEci(:,j) - offEci(:,centerGmm-(conj_n-1)*ord);
        for i = 1:n
            eta2 = pp.e2b(2,:,kk)';
            d    = dot(ell(:,i,kk),eta1)/dot(eta1,eta2)*eta2;
            p    = e2b*(ell(:,i,kk) - d - xc);
            X    = p(1)*cos(alpha(j)) + p(2)*sin(alpha(j));
            Y    = p(1)*sin(alpha(j)) - p(2)*cos(alpha(j));
            checkEll = sum([X;Y].^2./a(:,j).^2);
            if checkEll <= 1
                ell(:,i,kk) = nan;
            end
        end
    end
end
ellNew = squeeze(reshape(ell,[3,n*ord,1]));

% find closest point to original
for j = 1:ord
    pOrig     = rEci-center;
    [~,b]     = min(normOfVec(ellNew - pOrig));
    jj        = j + (conj_n - 1)*ord;
    distEta   = pp.e2b(2,:,jj)*pOrig;
    aa        = ellNew(:,b) + distEta*pp.e2b(2,:,jj)';
    pointEci(:,j) = aa + center;
end
return
%% Debug plot
plot3(ell(1,:,1),ell(2,:,1),ell(3,:,1),'LineWidth',2)
hold on
plot3(ell(1,:,2),ell(2,:,2),ell(3,:,2),'LineWidth',2)
plot3(ell(1,:,3),ell(2,:,3),ell(3,:,3),'LineWidth',2)
% plot3(ell(1,:,4),ell(2,:,4),ell(3,:,4),'LineWidth',2)
% plot3(ell(1,:,5),ell(2,:,5),ell(3,:,5),'LineWidth',2)
ppp = rEci-center;
plot3(ppp(1),ppp(2),ppp(3),'*')
plot3(aa(1),aa(2),aa(3),'*')
%%
legend('mixand 1','mixand 2','mixand 3','Original point','Solution')
grid on
xlabel('\xi')
ylabel('\zeta')
xlabel('\xi [-]')
ylabel('\zeta [-]')
end