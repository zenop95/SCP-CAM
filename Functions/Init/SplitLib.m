function [mu_i,w_i] = SplitLib(N)

load('uni_mat.mat','mu','w') % mu, w

mu_i    = mu(string(N));
w_i     = w(string(N));