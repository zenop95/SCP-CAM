function y = iteratePages(fun,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N = size(x,3);
y = nan(size(x));
for i = 1:N
    y(:,:,i) = fun(x(:,:,i));
end
end