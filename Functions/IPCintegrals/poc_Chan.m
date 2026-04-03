function Poc = poc_Chan(R,P,smd,varargin)
% Compute the collision probability between two objects using
% Chan's Method (REF. "Spacecraft Collision Probability", F. Kenneth Chan)
% -------------------------------------------------------------------------
% INPUTS: 
% - INPUT  : struct with the following fields: [m], [rad], [s]
%           (a_1,e_1,th_1,Rc,T,sa,re,sigma_ksi,sigma_zeta,rho,chi,phi,psi,u)
% - delta_r: relative position in the Bplane
% -------------------------------------------------------------------------
% OUTPUT: 
% - Poc : Collision Probability
% -------------------------------------------------------------------------
% Author:   Maria Francesca Palermo, Politecnico di Milano, 27 November 2020
%           e-mail: mariafrancesca.palermo@mail.polimi.it
if nargin > 3; ord = varargin{1}; else; ord = 5; end
sigma_xi   = sqrt(P(1,1));
sigma_zeta = sqrt(P(2,2));
rho        = abs(P(1,2))/(sigma_xi*sigma_zeta);
u          = R^2/(sigma_xi*sigma_zeta*sqrt(1.0-rho^2));
v          = smd;


% -------------------------------------------------------------------------
% The series can be truncated for m = 3 for small values of u.

m_vect = 0:ord; % hp) small values of u
P_mpart = 0;

for j = 1:length(m_vect)
    m = m_vect(j);
    k_vect = 0 : 1 : m_vect(j);
    P_kpart = 0;
    for w = 1:length(k_vect)
        k = k_vect(w);
        P_kpart = P_kpart + u^k/(2^k*factorial(k));
    end
    P_mpart =  P_mpart + (  v^m/(2^m* factorial(m))* (1- exp(-u/2)*P_kpart ));
end

% Probability of collision
Poc = exp(-v/2) * P_mpart;

end

