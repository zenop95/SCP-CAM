function q = rmat2quat(R)
% Compute the quaternion corresponding to the input rotation matrix
%
% Usage:
%   q = rmat2quat(R)
%
% Inputs:
%   R: 3x3xN rotation matrix
% Outputs:
%   q: 4xN corresponding  unit quaternion

validateattributes(R,{'numeric'},{'size',[3,3,NaN],'real','finite','nonnan'})

n = size(R, 3);
% Switch approach depending on input size
% Required by Simulink Code Generation
if n > 1
    q = cell2mat(arrayfun(@(k) rmat2quat1D(R(:,:,k)), 1:n, 'UniformOutput', false));
else
    q = rmat2quat1D(R);
end
end

function q = rmat2quat1D(R)
q    = zeros(4,1);
threshold = 1e-12;

tr = R(1,1) + R(2,2) + R(3,3);
if (tr > threshold) % if tr > 0:
    % scalar part:
    q(1,1) = 0.5*sqrt(tr + 1);
    s_inv = 1/(q(1,1)*4);
    % vector part:
    q(2,1) = (R(3,2) - R(2,3))*s_inv;
    q(3,1) = (R(1,3) - R(3,1))*s_inv;
    q(4,1) = (R(2,1) - R(1,2))*s_inv;
else % if tr <= 0, find the greatest diagonal element for calculating
    % the scale factor s and the vector part of the quaternion:
    if ( (R(1,1) > R(2,2)) && (R(1,1) > R(3,3)) )
        q(2,1) = 0.5*sqrt(R(1,1) - R(2,2) - R(3,3) + 1);
        s_inv = 1/(q(2,1)*4);
        
        q(1,1) = (R(3,2) - R(2,3))*s_inv;
        q(3,1) = (R(2,1) + R(1,2))*s_inv;
        q(4,1) = (R(3,1) + R(1,3))*s_inv;
    elseif (R(2,2) > R(3,3))
        q(3,1) = 0.5*sqrt(R(2,2) - R(3,3) - R(1,1) + 1);
        s_inv = 1/(q(3,1)*4);
        
        q(1,1) = (R(1,3) - R(3,1))*s_inv;
        q(2,1) = (R(2,1) + R(1,2))*s_inv;
        q(4,1) = (R(3,2) + R(2,3))*s_inv;
    else
        q(4,1) = 0.5*sqrt(R(3,3) - R(1,1) - R(2,2) + 1);
        s_inv = 1/(q(4,1)*4);
        
        q(1,1) = (R(2,1) - R(1,2))*s_inv;
        q(2,1) = (R(3,1) + R(1,3))*s_inv;
        q(3,1) = (R(3,2) + R(2,3))*s_inv;
    end
end
end