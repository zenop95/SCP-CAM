function [ln_q] = quat_log(q)
%Computes the logarithm of a quaternion
% Input:
%   q:       4-by-M matrix containing M quaternions.
%
% Output:
%   e_q:     4-by-M matrix containing natural logarithm of quaternions.

validateattributes(q,{'numeric'},{'nrows',4,'real','finite','nonnan'})
q = quat_normalize(q);

len = size(q,2);
v_norm = arrayfun(@(i) norm(q(2:4,i)),1:len, 'UniformOutput',true);

ln_q = reshape(cell2mat(arrayfun(@(i)...
    [0; log(norm(q(:,i)))+q(2:4,i)/v_norm(i)*acos(q(1,i)/norm(q(:,i)))],...
    1:len, 'UniformOutput',false)'),4,len);

ln_q(isnan(ln_q)) = 0;
end

