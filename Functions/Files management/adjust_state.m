function [mean_state_ECI, covariance_ECI] = adjust_state(mean_state_input,covariance_RTN,ref_frame,TCA_et)


if strcmp(ref_frame,'ITRF')
    
    % Rotation matrix from ECEF to ECI
    R_ecef2eci=cspice_sxform('ITRF93','J2000',TCA_et);
    
    % Mean state rotation from ECEF to ECI
    mean_state_ECI = R_ecef2eci*mean_state_input;
else
    mean_state_ECI = mean_state_input;
end

% Covariance rotation
Rot_eci2rtn = build_R_eci2rtn(mean_state_ECI);
Rot_eci2rtn = Rot_eci2rtn.';
R_rtn2eci_ext = blkdiag(Rot_eci2rtn,Rot_eci2rtn);

covariance_ECI=R_rtn2eci_ext*covariance_RTN*R_rtn2eci_ext.';


end

