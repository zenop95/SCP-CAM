function quat_res = quat_inner_product(q_1,q_2)
% Compute the inner product of the two input quaternions
% q_1 = [eta_1; eps_1], q_2 = [eta_2; eps_2] -> quat_res = [eta_1*eta_2+esp_1'*eps_2; zeros(3,1)]
%
%   Input:
%   q_1:            4-by-M matrix containing M quaternions.
%   q_2:            4-by-M matrix containing M quaternions.

%   Output:
%   quat_res:     4-by-M matrix containing inner product of quaternions of quaternions.

validateattributes(q_1,{'numeric'},{'nrows',4,'real','finite','nonnan'})
validateattributes(q_2,{'numeric'},{'nrows',4,'real','finite','nonnan'})

len1 = size(q_1,2); 
len2 = size(q_2,2);

if len1 ~= len2
    error('The dimensions of the array of q_1 and those of q_2 are not consistent')
end

quat_res = reshape(cell2mat(arrayfun(@(i) [q_1(1,i)*q_2(1,i)+q_1(2:4,i)'*q_2(2:4,i);...
    zeros(3,1)],1:len1, 'UniformOutput',false)'),4,len1);