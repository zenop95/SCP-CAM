function eulang = rmat2eulang(R, seq)
% Compute the the Euler angles defining the input rotation matrix according
% to the provided sequence (e.g.,'xyz', 'zyx')
%
% Usage:
%   eulang = rmat2eulang(R, seq)
%
% Inputs:
%   R: 3x3xN rotation matrix
%   seq: (string) angles sequence
% Outputs:
%   eul_ang: 3xN vector defining the corresponding Euler angles [rad]
%   Euler's angles are returned in the order specified by the input sequence

validateattributes(R,{'numeric'},{'size',[3,3,NaN],'real','finite','nonnan'})
validateattributes(seq,{'char'},{'nonempty'})

s = lower(seq) - 'w';
assert(all(s < 4 & s > 0), 'Invalid sequence ''%s''.', seq)

n = size(R, 3);
% Switch approach depending on input size
% Required by Simulink Code Generation
if n > 1
    eulang = cell2mat(arrayfun(@(k) rmat2eulang1D(R(:,:,k), s), 1:n,...
        'UniformOutput', false));
else
    eulang = rmat2eulang1D(R, s);
end
end

function eulang = rmat2eulang1D(R, s)
eulang = zeros(3,1);
R = R';
g = 2*mod(s(1)-s(2),3)-3;   % +1 -> seq(2)>seq(1),  -1 -> seq(1)>seq(2)
h = 2*mod(s(2)-s(3),3)-3;   % +1 -> seq(3)>seq(2),  -1 -> seq(2)>seq(3)
[sin, cos, ~, ang] = atan2sinCos(h*R(s(2),s(1)),R(6-s(2)-s(3),s(1)));
R_tmp = R;
idx = 1+mod(s(3)+(0:1),3);
R_tmp(idx,:) = [cos sin; -sin cos]*R(idx,:);

eulang(1) = atan2(-g*R_tmp(s(2),6-s(1)-s(2)),R_tmp(s(2),s(2)));
eulang(2) = atan2(-g*R_tmp(6-s(1)-s(2),s(1)),R_tmp(s(1),s(1)));
eulang(3) = ang;
eulang = -eulang;

if ((s(1)~=s(3) && abs(eulang(2))>pi/2)) % || seq(1)==seq(3) && eul_ang(2)<0))  % remove redundancy
    mk          =   s(1)~=s(3);
    eulang(2)  =   (2*mk-1)*eulang(2);
    eulang     =   eulang-((2*(eulang>0)-1) .* [1; mk; 1])*pi;
end
end