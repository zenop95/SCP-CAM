function qn = quat_normalize(q)
% Compute the normalized  quaternion
% Input:
%   q:  4xN matrix containing M quaternions.
% Output:
%   qn: 4xN matrix containing M normalized quaternions.

len = size(q,2);
validateattributes(q,{'numeric'},{'nrows',4,'real','finite','nonnan'})

% Normalize every quaternion and concatenate over 2nd dimension
qn = reshape(cell2mat(arrayfun(@(i) q(:,i)./quat_norm(q(:,i)),...
    1:len, 'UniformOutput', false)), 4, len);
end
