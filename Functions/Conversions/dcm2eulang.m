function eulang = dcm2eulang(dcm, seq)
% Compute the the Euler angles defining the input direction cosine matrix 
% according to the provided sequence (e.g.,'xyz', 'zyx')
%
% Usage:
%   eulang = dcm2eulang(dcm, seq)
%
% Inputs:
%   R: 3x3xN direction cosine matrix. The provided DCM is in the
%      pre-multiply format (i.e., v2 = dcm21 * v1).
%   seq: (string) angles sequence
%
% Outputs:
%   eul_ang: 3xN vector defining the corresponding Euler angles [rad]
%   Euler's angles are returned in the order specified by the input sequence

validateattributes(dcm,{'numeric'},{'size',[3,3,NaN],'real','finite','nonnan'})
validateattributes(seq,{'char'},{'nonempty'})

s = lower(seq) - 'w';
assert(all(s < 4 & s > 0), 'Invalid sequence ''%s''.', seq)

% Transpose each dcm to obtain rotation matrices
R = permute(dcm, [2 1 3]);

% Compute Euler's angles
eulang = rmat2eulang(R, seq);
if any(isnan(eulang)); warning('nan value in eulang set to 0'); end
eulang(isnan(eulang)) = 0;
end