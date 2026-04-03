function [] = plotEllipsoidPoints(ellipsoids,ellipsoids2)
%plotEllipsoidPoints Plots the covariance ellipsoid for the required node
%and the interior point or the successive convexification points.

point   = ellipsoids.z;
Dr      = ellipsoids.Dr;
cov2CW  = ellipsoids.cov2CW;
a       = ellipsoids.semiaxes(1);
b       = ellipsoids.semiaxes(2);
c       = ellipsoids.semiaxes(3);
DrNew   = ellipsoids2.Dr;

% Find the plane tangent to the ellipsoid on the indentified point
gradNorm = ellipsoids.grad/norm(ellipsoids.grad);
[XX, YY] = meshgrid(-b:b/10:b); 
ZZ = -1/gradNorm(3)*(gradNorm(1)*XX + gradNorm(2)*YY);
tanPlane(:,:,1) = XX + point(1);
tanPlane(:,:,2) = YY + point(2);
tanPlane(:,:,3) = ZZ + point(3);

[ax,ang] = quat2axang(dcm2quat(cov2CW));
[X,Y,Z] = ellipsoid(0,0,0,a,b,c,1000);
S = surf(X,Y,Z,'FaceAlpha',0.1);
shading interp
rotate(S,ax,-ang*180/pi)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
hold on
surf(tanPlane(:,:,1),tanPlane(:,:,2),tanPlane(:,:,3),'FaceAlpha',0.3);
shading interp
plot3([0 Dr(1)],[0 Dr(2)],[0 Dr(3)],'-o')
plot3([0 point(1)],[0 point(2)],[0 point(3)],'-*')
plot3([0 DrNew(1)],[0 DrNew(2)],[0 DrNew(3)],'-o')
legend('','','Original point','Point on Ellipsoid','New Point')
axis equal
hold off
end