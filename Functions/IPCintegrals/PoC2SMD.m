function [SMD, err, iter] = PoC2SMD(P, R, PoC_enforced, m, n,  varargin)
% Analytic SMD computation given an enforced PoC with M =1
% -------------------------------------------------------------------------
% INPUTS:
% - P = summed Covariance matrix of the primary and the secondary
% - R = summed primary and secondary object radii
% - PoC_enforced : the desired collision probability
% - m : order accuracy
% - n : order tayor polynomial to approximate exponential
% - varargin{1}  : pecerntual error wrt the previous SMD estimation. This
%                  variable enables an iterative process until the required
%                  error on the SMD is met.
% - varargin{2}  : maximum number of desired iterations. if varargin{2} is
%                  not set, then this variable is 20. if
%                  varargin{1} is empty as well, then the number of
%                  iteration is 1.
% -------------------------------------------------------------------------
% OUTPUT:
% - SMD = Output Square Mahalanobis Distance
%
% -------------------------------------------------------------------------
% Author: Andrea De Vittori, Politecnico di Milano, 04 March 2022

% optional variables management
len_varargin = length(varargin);
if len_varargin >= 1
    err_threshold = varargin{1};
else
    err_threshold = inf;
end

if len_varargin == 2
    max_iter = varargin{2};
else
    max_iter = 20;
end
% Equivalent area
u = R^2/(sqrt(P(1, 1)*P(2, 2))*sqrt(1-P(1, 2)^2/(P(1, 1)*P(2, 2))));

% SMD at order 0
v0 = -2*log(PoC_enforced/(1-exp(-u/2)));

% variables initialization
iter = 0;
err = inf;
if m == 0
    SMD= v0;
    err = 0;
elseif  m >= 1
    dim_taylor = n + 1;
    dim_chan = m + 1;
    dim_tot = m + n +1;
    err_vect = ones(2,1)*inf;
    c_i = zeros(1, dim_tot);
    b_i = zeros(1, dim_chan);
    a_i = zeros(1, dim_taylor);
    % Chan's PoC at any order
    while err >= err_threshold && iter < max_iter

        % Build the Taylor coefficients
        coeff_exp = exp(-v0/2);
        err_vect(1) = err_vect(2);
        for i=0:n
            a_n = 0;
            for p = i:n
                 a_n  =  a_n +(-0.5)^p/factorial(p)*nchoosek(p, i)*(-v0)^(p-i);
            end
            a_i(i+1)=   coeff_exp*a_n ;
        end

        % Build the Chan's coefficients
        s_k = zeros(2,1);
        S_k = zeros(2,1);
        s_k(1) = exp(-u/2);
        S_k(1) = s_k(1);
        b_i(1) = 1 - S_k(1);
        for k = 1:m
            s_k(2) = u/(2*k)*s_k(1);
            S_k(2) = S_k(1) + s_k(2);
            q_m = 1 - S_k(2);
            b_i(k+1) = 1/(2^k*factorial(k))*q_m;
            s_k(1) = s_k(2);
            S_k(1) = S_k(2);
        end

        for j=0:n+m
            c_nm =0;

            if n>m
                ind_sup = min(n, j);
                ind_inf = max(0, j-m);
                for i =ind_inf:ind_sup

                    c_nm =  c_nm+b_i(j-i+1)*a_i(i+1);


                end
            else
                ind_sup = min(m, j);
                ind_inf = max(0, j-n);
                for i =ind_inf:ind_sup
                    c_nm =  c_nm+b_i(i+1)*a_i(j-i+1);
                end
            end
            c_i(j+1) = c_nm;
        end
        c_i(1) = c_i(1) - PoC_enforced;
        SMD = roots(flip(c_i));
        SMD = SMD(SMD == real(SMD)& SMD >0 & abs(v0-SMD)== min(abs(v0-SMD)));
        if isempty(SMD)
            SMD = v0;
            err = err_vect(1);
            break
        end
        err_vect(2) = abs((SMD-v0)/v0)*100;
        err = err_vect(2);
        v0 = SMD;
        iter = iter + 1;
    end
end

end




