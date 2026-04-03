function [conjunction] = objects_definition(CDM_input_struct)

conjunction.mu = 398600.4415; % Earth gravity constant [km^3/s^2]


% ------------------------------- TCA epoch -------------------------------

% TCA in ephemeral time
refdate_TCA = CDM_input_struct.TCA;
conjunction.TCA_et = cspice_str2et(char(refdate_TCA));

% ------------------------------- Object 1 --------------------------------
object1_mean_state = CDM_input_struct.object1.data.mean_state;
R_p = CDM_input_struct.object1.data.radius;
% Covariance in RTN
object1_covariance_RTN = CDM_input_struct.object1.data.covariance;
% Mean state reference frame
object1_RF = CDM_input_struct.object1.reference_frame;

[object1_mean_state_ECI_TCA, object1_covariance_ECI] = adjust_state(...
    object1_mean_state,object1_covariance_RTN,object1_RF,conjunction.TCA_et);

% Save object1 in a dedicated structure

conjunction.object1 = CDM_input_struct.object1;
conjunction.object1.mean_state = object1_mean_state_ECI_TCA;
conjunction.object1.covariance_pos = object1_covariance_ECI(1:3,1:3);
conjunction.object1.radius = R_p;

% ------------------------------- Object 2 --------------------------------
object2_mean_state = CDM_input_struct.object2.data.mean_state;
R_s = CDM_input_struct.object2.data.radius;
% Covariance in RTN
object2_covariance_RTN = CDM_input_struct.object2.data.covariance;
% Mean state reference frame
object2_RF = CDM_input_struct.object2.reference_frame;

[object2_mean_state_ECI_TCA, object2_covariance_ECI] = adjust_state(...
    object2_mean_state,object2_covariance_RTN,object2_RF,conjunction.TCA_et);

% Save object2 in a dedicated structure
conjunction.object2 = CDM_input_struct.object2;
conjunction.object2.mean_state = object2_mean_state_ECI_TCA;
conjunction.object2.covariance_pos = object2_covariance_ECI(1:3,1:3);
conjunction.object2.radius = R_s;

% -------------------------------------------------------------------------

R_earth = 6371.01; % Earth radius [km]
g0 =conjunction.mu/R_earth^2; % Earth gravity acceleration [km/s^2]

% Primary propulsion parameters (IRIDIUM)
conjunction.object1.Tmax = 0.1; % [N] thrust
conjunction.object1.Isp = 220; % [s]
conjunction.object1.m0 = 500; % [kg] mass
conjunction.object1.ce = conjunction.object1.Isp*g0;% effective velocity [km/s]

end

