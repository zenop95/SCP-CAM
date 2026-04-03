function [] = plotEllipsoidInit(ellipsoids)
%plotEllipsoidInit Plots the covariance ellipsoid for the required node
%and the interior point or the successive convexification points.

point   = ellipsoids.z;
Dr      = ellipsoids.Dr;
cov2CW  = ellipsoids.cov2CW;
a       = ellipsoids.a;
b       = ellipsoids.b;
c       = ellipsoids.c;

[ax,ang] = quat2axang(dcm2quat(cov2CW));
[X,Y,Z] = ellipsoid(0,0,0,a,b,c,1000);
S = surf(X,Y,Z,'FaceAlpha',0.1);
shading interp
rotate(S,ax,-ang*180/pi)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
hold on
plot3(Dr(1),Dr(2),Dr(3),'o')
plot3(point(1),point(2),point(3),'o')
legend('','Original point','Point on Ellipsoid')
axis equal
hold off
end