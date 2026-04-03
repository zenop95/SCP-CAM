function eulang = quat2eulang(q, seq)
% Compute the the Euler angles defining the input rotation matrix according
% to the provided sequence (e.g.,'xyz','zyx')
%
% Usage:
%   eulang = quat2eulang(q, seq)
%
% Inputs
%   q: 4xN unit quaternion
%   seq: (string) angles sequence
%
% Outputs:
%   eulang: 3xN vector defining the corresponding Euler angles [rad]
%   Euler's angles are returned in the order specified by the input sequence

validateattributes(q,{'numeric'},{'nrows',4,'real','finite','nonnan'})
validateattributes(seq,{'char'},{'nonempty'})

s = lower(seq) - 'w';
assert(all(s < 4 & s > 0), 'Invalid sequence ''%s''.', seq)

eulang = rmat2eulang(quat2rmat(q), seq);
end