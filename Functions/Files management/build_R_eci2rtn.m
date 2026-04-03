function R_eci2rtn = build_R_eci2rtn(mean_state)
% Build up the rotation matrix from ECI to RTN reference frame
% -------------------------------------------------------------------------
% INPUTS:
% - mean_state: orbital mean state in ECI reference frame (6xn_rotations)
% -------------------------------------------------------------------------
% OUTPUT:
% - R_eci2rtn : matrix to rotate one or multiple vectors from ECI to RTN
%               reference frame (3x3xn_rotations)
%     - x-axis: from [1 0 0] to [position_unitary_vector 0 0]
%     - y-axis: from [0 1 0] to cross(z_axis,x_axis)
%     - z-axis: from [0 0 1] to [angular_momentum_unitary vector 0 0]
% -------------------------------------------------------------------------
% Author:   Marco Felice Montaruli, Politecnico di Milano, 24 August 2020
%           e-mail: marcofelice.montaruli@polimi.it


% Check that the mean state is correctly passed
if size(mean_state,1)~=6
    ME = MException('PoliMiSST:build_R_eci2rtn:inputError', ...
        'The input state has to be passed as a 6xn_rotations vector.');
    throw(ME)
end

% Number of rotations
n_rotations = size(mean_state,2);

% ECI position, velocity and angular momentum
position_eci = mean_state(1:3,:);
velocity_eci = mean_state(4:6,:);
angular_momentum_eci = cross(position_eci,velocity_eci);

% RTN unitary vectors
p_eci_vers = position_eci./vecnorm(position_eci);
h_eci_vers = angular_momentum_eci./vecnorm(angular_momentum_eci);
n_vect = cross(h_eci_vers,p_eci_vers);
n_vers = n_vect./vecnorm(n_vect);

% Rotation matrix from ECI to RTN
R_eci2rtn = zeros(3,3,n_rotations);
R_eci2rtn(1,:,:) = p_eci_vers;
R_eci2rtn(2,:,:) = n_vers;
R_eci2rtn(3,:,:) = h_eci_vers;


end

