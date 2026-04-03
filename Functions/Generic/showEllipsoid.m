function [] = showEllipsoid(c,P)
%plotEllipsoidPoints Plots the covariance ellipsoid for the required node
%and the interior point or the successive convexification points.
[V,D] = eig(P);
[a,b] = sort(diag(D),'descend');
D = diag(a);
cov2Eci = V(:,b);
p        = sqrt(diag(D));
[ax,ang] = quat2axang(dcm2quat(cov2Eci));
[X,Y,Z] = ellipsoid(c(1),c(2),c(3),p(1),p(2),p(3),100);
S = surf(X,Y,Z,'FaceAlpha',0.1);
shading interp
% rotate(S,ax,-ang*180/pi)
hold on
plot3(c(1),c(2),c(3),'*')
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

end