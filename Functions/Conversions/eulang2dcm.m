function dcm = eulang2dcm(eulang, seq)
% Compute the direction cosine matrix corresponding to the input Euler 
% angles according to the provided sequence (e.g., 'xyz', 'zyx')
%
% Usage:
%   dcm = eulang2dcm(eulang, seq)
%
% Inputs:
%   eulang: 3xN vector defining the Euler angles according to the given sequence
%           Euler's angles must be input in the order specified by the
%           input sequence
%   seq: (string) angles sequence
%
% Outputs:
%   dcm: 3x3xN corresponding direction cosine matrix. The provided DCM is 
%       in the pre-multiply format (i.e., v2 = dcm21 * v1).

validateattributes(seq,{'char'},{'nonempty'})
validateattributes(eulang,{'numeric'},{'nrows',3,'real','finite','nonnan'})

s = lower(seq) - 'w';
assert(all(s < 4 & s > 0), 'Invalid sequence ''%s''.', seq)

R = eulang2rmat(eulang, seq);

% Transpose each rotmat to obtain a dcm
dcm = permute(R, [2 1 3]);
end