function [grad] = smd_grad(P,r)
% SmdGrad algebraically computes the gradient of the Squared Mahalanobis 
% Distance associated with the covariance matrix C and the value r of 
% the random variable. Obtained with the symbolic toolbox of MATLAB.
%
% INPUT: C = [-] (3,3) Covariance matrix
%        r = [-] (3,1) Relative position
%
% OUTPUT: grad = [-] (3,1) Gradient of the Squared Mahalanobis Distance
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
x = r(1); y = r(2); z = r(3);
p11 = P(1,1); p22 = P(2,2); p33 = P(3,3);
p12 = P(1,2); p13 = P(1,3); p23 = P(2,3);

t2 = p12.*p13;
t3 = p11.*p23;
t4 = p12.*p23;
t5 = p13.*p22;
t6 = p13.*p23;
t7 = p12.*p33;
t8 = p12.^2;
t9 = p13.^2;
t10 = p23.^2;
t11 = p11.*p22.*p33;
t12 = -t3;
t13 = -t5;
t14 = -t7;
t15 = p23.*t3;
t16 = p13.*t5;
t17 = p12.*t7;
t18 = p23.*t2.*2.0;
t20 = -t11;
t19 = -t18;
t21 = t2+t12;
t22 = t4+t13;
t23 = t6+t14;
t24 = t15+t16+t17+t19+t20;
t25 = 1.0./t24;
grad = [t25.*(t10.*x-t6.*y+t7.*y-t4.*z+t5.*z-p22.*p33.*x)+t25.*x.*(t10-p22.*p33)-t23.*t25.*y-t22.*t25.*z;-t25.*(t6.*x+t14.*x-t9.*y+t2.*z+t12.*z+p11.*p33.*y)+t25.*y.*(t9-p11.*p33)-t23.*t25.*x-t21.*t25.*z;-t25.*(t4.*x+t13.*x+t2.*y+t12.*y-t8.*z+p11.*p22.*z)+t25.*z.*(t8-p11.*p22)-t22.*t25.*x-t21.*t25.*y];