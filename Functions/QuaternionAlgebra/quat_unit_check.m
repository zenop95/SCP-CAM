function isUnit = quat_unit_check(q, varargin)
% Check wheter the input quaternion is unitary up to an imposed threshold
% quat = [eta; eps] -> +1 if -thresold <= eta^2+eps'*eps <= +threasold
%
% INPUT
%   q: 4xN matrix containing M quaternions.
%   threshold: scalar threshold (tollerance)
% OUTPUT
%   isUnit: 1xN boolean -> 1 = yes, 0 = no

validateattributes(q,{'numeric'},{'nrows',4,'real','finite','nonnan'})

threshold = 1e-8;
if nargin > 1; threshold = varargin{2}; end

len = size(q,2);
isUnit  = reshape(cell2mat(arrayfun(@(i) checkUnit(q(:,i),threshold),...
    1:len, 'UniformOutput', false)), 1, len);
end

function bool = checkUnit(q, threshold)
bool = (abs(quat_norm(q)-1) <= threshold);
end