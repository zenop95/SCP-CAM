function R = eulang2rmat(eulang, seq)
% Compute the rotation matrix corresponding to the input Euler angles
% according to the provided sequence (e.g., 'xyz', 'zyx')
%
% Usage:
%   R = eulang2rmat(eulang, seq)
%
% Inputs:
%   eulang: 3xN vector defining the Euler angles according to the given sequence
%           Euler's angles must be input in the order specified by the
%           input sequence
%   seq: (string) angles sequence
% Outputs:
%   R: 3x3xN corresponding rotation matrix

validateattributes(seq,{'char'},{'nonempty'})
validateattributes(eulang,{'numeric'},{'nrows',3,'real','finite','nonnan'})

s = lower(seq) - 'w';
assert(all(s < 4 & s > 0), 'Invalid sequence ''%s''.', seq)

n = size(eulang, 2);
% Switch approach depending on input size
% Required by Simulink Code Generation
if n > 1
    R = cell2mat(arrayfun(@(k) eulang2rmat1D(eulang(:,k), s), 1:n,...
        'UniformOutput', false));
else
    R = eulang2rmat1D(eulang, s);
end

% Concatenate over 3rd dimension
R = reshape(R, 3, 3, []);
end

function R = eulang2rmat1D(eulang, s)
E = eye(3);
R = axang2rmat(...
    E(:,s(1)),...
    eulang(1))*axang2rmat(E(:,s(2)),...
    eulang(2))*axang2rmat(E(:,s(3)),...
    eulang(3));
end