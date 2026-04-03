function q = eulang2quat(eulang, seq)
% Compute the rotation matrix corresponding to the input Euler angles
% according to the provided sequence (e.g., 'xyz', 'zyx')
%
% Usage:
%   q = eulang2quat(eulang, seq)
%
% Inputs:
%   eulang: 3xN vector defining the Euler angles according to the given sequence
%           Euler's angles must be input in the order specified by the 
%           input sequence
%   seq: (string) angles sequence
% Outputs:
%   q: 4xN corresponding unit quaternion

validateattributes(seq,{'char'},{'nonempty'})
validateattributes(eulang,{'numeric'},{'nrows',3,'real','finite','nonnan'})

s = seq - 'w';
assert(all(s < 4 & s > 0), 'Invalid sequence ''%s''.', seq)

R = eulang2rmat(eulang, seq);
q = rmat2quat(R);
end