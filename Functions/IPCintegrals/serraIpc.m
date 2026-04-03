function IPC = serraIpc(mu,P,R,n)
% cuboidIpc computes the Instantaneous Collision Probability (IPC) as the
% integral of the gaussian probability distribution ellipsoid of the 
% relative position of two bodies described by the mean and covariance 
% matrix in input over the Hard Body Sphere of radius R

% INPUT: mu = [m] expected value of the relative position random variable
%        P  = [m^2] Covariance Matrix of the relative position random variable
%        R  = [m] Hard Body Radius (HBR) of the collision

% OUTPUT: Ipc = [-] Computed Instantaneous Collision Probability


% Bibliography: Zhang, S., Fu, T., Chen, D., & Cao, H. (2020). 
%               Satellite instantaneous  collision probability computation 
%               using equivalent volume cuboids. Journal of Guidance, 
%               Control, and Dynamics, 43(9), 1757–1763. 

% Author: Zeno Pavanello, 2022
%--------------------------------------------------------------------------

%Factorize covariance matrix
%Find orthogonal matrix of eigenvectors of A and its eigenvalues
[V,D] = eig(P);
[sigma,ord] = sort(sqrt(diag(D)),'ascend');
mu    = V(:,ord)'*mu;
sigma_1 = sigma(1); 
sigma_2 = sigma(2); 
sigma_3 = sigma(3); 
p       = 1/(2*sigma_1^2);
gamma_2 = p - 1/(2*sigma_2^2);
gamma_3 = p - 1/(2*sigma_3^2);
gamma_a = gamma_2 + gamma_3;
gamma_m = gamma_2*gamma_3;
theta_1 = mu(1)^2/(2*sigma_1^4);
theta_2 = mu(2)^2/(2*sigma_2^4);
theta_3 = mu(3)^2/(2*sigma_3^4);

q_0     = -2;
q_1     = 4*gamma_a + 2*p;
q_2     = -(2*gamma_a^2 + 4*p*gamma_a + 4*gamma_m);
q_3     = (2*gamma_a^2 + 4*gamma_m)*p + 4*gamma_a*gamma_m;
q_4     = -4*p*gamma_a*gamma_m - 2*gamma_m^2;
q_5     = 2*p*gamma_m^2;
q_6     = 0;
q       = [q_0; q_1; q_2; q_3; q_4; q_5; q_6];
f_0     =  -theta_3 - theta_2 - gamma_a;
f_1     = 2*theta_3*gamma_2 + 2*theta_2*gamma_3 + gamma_a^2 + 2*gamma_m;
f_2     = -theta_3*gamma_2^2 - theta_2*gamma_3^2 - 3*gamma_a*gamma_m;
f_3     = 2*gamma_m^2;
s_0     = f_0 - 2*p;
s_1     = (-f_0 + 4*gamma_a)*p + f_1;
s_2     = (-f_1 - 4*gamma_m - 2*gamma_a^2)*p + f_2;
s_3     = (-f_2 + 4*gamma_m*gamma_a)*p + f_3;
s_4     = -4*gamma_m^2*p;
s_5     = 0;
s       = [s_0; s_1; s_2; s_3; s_4; s_5];
omega   = nan(6,1);
for i = 1:6
    omega(i) = theta_1/2*q(i)+s(i);
end
C = exp(-dot(mu.^2,1./(2*sigma.^2)))/(2^(3/2)*prod(sigma));
alpha_0 = C;
alpha_1 = -omega(1)*alpha_0/2;
alpha_2 = ((q_1-omega(1))*alpha_1 - omega(2)*alpha_0)/4;
alpha_3 = ((2*q_1-omega(1))*alpha_2 + ...
             (q_2-omega(2))*alpha_1 - omega(3)*alpha_0)/6;
alpha_4 = ((3*q_1-omega(1))*alpha_3 + (2*q_2-omega(2))*alpha_2 + ...
             (q_3-omega(3))*alpha_1 - omega(4)*alpha_0)/8;
alpha_5 = ((4*q-1-omega(1))*alpha_4 + (3*q_2 - omega(2))*alpha_3 + ...
           (2*q_3-omega(3))*alpha_2 + (q_4 - omega(4))*alpha_1 - ...
                  omega(5)*alpha_0)/10;
alpha   = [alpha_0; alpha_1; alpha_2; alpha_3; alpha_4; alpha_5];
c_0     = 4*alpha_0*R^3/(3*sqrt(pi));
c = nan(n,1);
c(1) = c_0;
for i = 1:5
    prodj = 1;
    for j = 0:i-1
        prodj = prodj*(j + 5/2);
    end
    c(i+1) = 4*alpha(i+1)*R^(2*i+3)/(3*sqrt(pi)*prodj);
end
for k = 7:n
    w = 0;
    for i = 1:6
        den = 1;
        for j = 1:i
            den = den*(k - 1 + 5/2 - j);
        end
        w = w + R^(2*i)/den*((k - 1 - i)*q(i+1) - omega(i))*c(k-i);
    end
    c(k) = 1/(2*(k - 1))*w;
end
IPC = exp(-p*R^2)*sum(c); 
IPC(IPC<0) = 0;
IPC(IPC>1) = 1;
end
