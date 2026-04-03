function smd = smd_sym(C,r)
% Smd algebraically computes the Squared Mahalanobis Distance associated 
% with the covariance matrix C and the value r of the random variable.
%
% INPUT: C = 3x3 covariance matrix
%        r = 3x1 position random variable
%
% OUTPUT: smd = Squared Mahalanobis Distance
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
x = r(1); y = r(2); z = r(3);
c11 = C(1,1); c22 = C(2,2); c33 = C(3,3);
c12 = C(1,2); c13 = C(1,3); c23 = C(2,3);

t2 = c12.^2;
t3 = c13.^2;
t4 = c23.^2;
t5 = c11.*c22.*c33;
t9 = c12.*c13.*c23.*2.0;
t6 = c11.*t4;
t7 = c22.*t3;
t8 = c33.*t2;
t10 = -t9;
t11 = -t5;
t12 = t6+t7+t8+t10+t11;
t13 = 1.0./t12;
smd = t13.*x.*(t4.*x-c22.*c33.*x-c13.*c23.*y+c12.*c33.*y-c12.*c23.*z+c13.*c22.*z)+t13.*y.*(t3.*y-c13.*c23.*x+c12.*c33.*x-c11.*c33.*y-c12.*c13.*z+c11.*c23.*z)+t13.*z.*(t2.*z-c12.*c23.*x+c13.*c22.*x-c12.*c13.*y+c11.*c23.*y-c11.*c22.*z);
