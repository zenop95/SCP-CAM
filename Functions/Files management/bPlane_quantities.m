function  [OUTPUT, R, K] = bPlane_quantities(conjunction,Kp,Ks, Chan_order)

% Compute the b-plane quantities and the kinematics matrix using 
% -------------------------------------------------------------------------
% INPUTS:
% - object1: structure for object1, whose fields are
%               - mean state at TCA [km, km/2]
%               - position covariance at TCA [km^2]
%               - radius [km]
% - object2: structure for object1, whose fields are
%               - mean state at TCA [km, km/2]
%               - position covariance at TCA [km^2]
%               - radius [km]
% - Kp    : Covariance Scale Factor for object 1
% - Ks    : Covariance Scale Factor for object 2
% - Rb : rotation matrix from ECI to Bplane
% - TCA_et : TCA in ephemeral time
% -------------------------------------------------------------------------
% OUTPUT:
% - OUTPUT   : struct with the following fields: [km]
%             (a_1,e_1,th_1,Rc,T,sa,re,C,sigma_ksi,sigma_zeta,rho,chi,phi,psi,u)
% -------------------------------------------------------------------------
% Author:   Maria Francesca Palermo, Politecnico di Milano, 27 November 2020
%           e-mail: mariafrancesca.palermo@mail.polimi.it




% -------------------------------------------------------------------------
% CONSTANTS
mu = 398600.4415; % Earth gravity constant [km^3/s^2]


% -------------------------------------------------------------------------
% INPUT MANAGEMENT (CDM -> INPUT)

% Objects' covariances in ECI r.f. @t_CA 
Cov1_eci = conjunction.object1.covariance_pos; % [km^2]
Cov2_eci = conjunction.object2.covariance_pos; % [km^2]

% Objects' position and velocity in ECI r.f. @t_CA

r_1_eci = conjunction.object1.mean_state(1:3); % [km]
v_1_eci = conjunction.object1.mean_state(4:6); % [km/s]

r_2_eci = conjunction.object2.mean_state(1:3); % [km]
v_2_eci = conjunction.object2.mean_state(4:6); % [km/s]


% Primary object orbital elements [km, rad]
elems_1 = car2elements(conjunction.object1.mean_state, mu);
elems_2 = car2elements(conjunction.object2.mean_state, mu);

OUTPUT.th_1 = elems_1(6); % True anomaly at TCA
OUTPUT.e_1 = elems_1(2); % Eccentricity
OUTPUT.a_1 = elems_1(1); % semi-major axis [km]

% Position norm of the first object @t_CA
OUTPUT.Rc = norm(r_1_eci); %[km]

% Orbital period of the first object 
OUTPUT.T = elems_1(7); % Orbital period [s]
OUTPUT.T_2 = elems_2(7); % Orbital period [s]

% Orbital angular momentum
h_1= cross(r_1_eci, v_1_eci);
uh_1= h_1/norm(h_1);

% Sum of the two objcets radii

if isfield(conjunction.object1,'radius') && isfield(conjunction.object2,'radius')
    OUTPUT.sa = conjunction.object1.radius+conjunction.object2.radius; %[km]
elseif isfield(conjunction,'sa')
    OUTPUT.sa = conjunction.sa;
end


% Nominal relative distance between the two objects @t_CA
re_eci = r_1_eci - r_2_eci; % [km] ECI r.f.

% Rotation matrix from ECI to Bplane
Rb = R_eci2bplane(conjunction.object1.mean_state(4:6),conjunction.object2.mean_state(4:6));
OUTPUT.Rb = Rb;

Rb_2D = [Rb(1,1) Rb(1,2)  Rb(1,3);
         Rb(3,1) Rb(3,2)  Rb(3,3)];
     
OUTPUT.re = Rb*re_eci; % [km]  B-Plane r.f.
OUTPUT.be = Rb_2D*re_eci; % [km]  B-Plane r.f.
OUTPUT.Rb_2D = Rb_2D;

% -------------------------------------------------------------------------
% Overall covariance projected in b-plane
Cov1_bp = Rb*Cov1_eci*Rb'; %[km^2]
Cov2_bp = Rb*Cov2_eci*Rb'; %[km^2]
Cov_bp = Cov1_bp.*Kp+Cov2_bp.*Ks; %[km^2]

% ksi-zeta submatrix
C = [Cov_bp(1,1) Cov_bp(1,3); 
    Cov_bp(3,1) Cov_bp(3,3)]; %[km^2]

OUTPUT.C = C;
OUTPUT.sigma_ksi = sqrt(C(1,1)); % std deviation along ksi [km]
OUTPUT.sigma_zeta =  sqrt(C(2,2)); % std deviation along zeta [km]
OUTPUT.sigmaa = sqrt(OUTPUT.sigma_ksi*OUTPUT.sigma_zeta); % [km]
OUTPUT.rho = C(1,2)/(OUTPUT.sigma_ksi*OUTPUT.sigma_zeta); % correlation factor [-]

% -------------------------------------------------------------------------
% Angles in b-plane [rad]
OUTPUT.chi = norm(v_2_eci)/norm(v_1_eci);
% Let the velocity of S2 at collision be related to the velocity of S1 by a (positive) rotation of angle
OUTPUT.phi =  atan2(dot(cross(v_1_eci,v_2_eci),uh_1), dot(v_1_eci,v_2_eci)); 
% followed by an out-of-plane rotation −π∕2 < ψ < π∕2 in the direction approaching uh1
OUTPUT.psi = atan( ((dot(v_2_eci,uh_1))*norm(cross(v_2_eci,uh_1))) / (dot(v_2_eci,v_2_eci) - dot(v_2_eci,uh_1)^2) ); 

% -------------------------------------------------------------------------
% ratio of the impact cross-sectional area to the area of the 1σ covariance ellipse in the b-plane
OUTPUT.u = OUTPUT.sa^2 / (OUTPUT.sigma_ksi*OUTPUT.sigma_zeta * sqrt(1-OUTPUT.rho^2)); 


% -------------------------------------------------------------------------
% Useful matrix
OUTPUT.Lx = [1 0 0; 0 0 0; 0 0 0];
OUTPUT.Lz = [0 0 0; 0 0 0; 0 0 1];
OUTPUT.S = [1 0 0 ; 0 0 0; 0 0 1];

% -------------------------------------------------------------------------

% Compute rotation and kinematics matrices R and K
[OUTPUT.R, OUTPUT.K] = kinematics(OUTPUT);

% rotation matrix from ECI to bplane
OUTPUT.R_2D = [OUTPUT.R(1,1:3); OUTPUT.R(3,1:3)];

% Chan order expansion
OUTPUT.Chan_order = Chan_order;

end

