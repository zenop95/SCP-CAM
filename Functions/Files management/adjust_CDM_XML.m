function conjunction_struct=adjust_CDM_XML(CDM_filename)
% Read the whole CDM content and adjust conjunction most remarkable
% parameters
% -------------------------------------------------------------------------
% INPUTS:
% - filename: CDM filename
% -------------------------------------------------------------------------
% OUTPUT:
% - conjunction_struct: strucuture containing conjunction most remarkable
%                       parameters. Structure fields are:
%                   - message_ID
%                   - TCA
%                   - miss_distance
%                   - object1_struct
%                       - reference frame
%                       - maneuverable (optional - if present)
%                       - metadata (optional - if present)
%                       - data.radius
%                       - data.mean_state
%                       - data.covariance
%                   - object2_struct
%                       - reference frame
%                       - maneuverable (optional - if present)
%                       - metadata (optional - if present)
%                       - data.radius
%                       - data.mean_state
%                       - data.covariance
% -------------------------------------------------------------------------
% Author:        Marco Felice Montaruli, Politecnico di Milano, 27 November 2020
%                e-mail: marcofelice.montaruli@polimi.it


% Read CDM file
CDM_struct=read_CDM_XML(CDM_filename);

% ----------------------- CDM general information -------------------------

% Message ID
conjunction_struct.message_ID = CDM_struct.header.contents.MESSAGE_ID;

% TCA
conjunction_struct.TCA = CDM_struct.body.relativeMetadataData.contents.TCA;

% Miss distance
conjunction_struct.miss_distance = str2double(CDM_struct.body.relativeMetadataData.contents.MISS_DISTANCE)*1e-3; % [km]


% -------------------------------- Object 1 -------------------------------

% Object 1 structure
object1_struct = build_obj_struct(CDM_struct.body.object1);

% Reference frame
object1_struct.reference_frame = CDM_struct.body.object1.metadata.contents.REF_FRAME;

% Save in the Conjunction structure
conjunction_struct.object1 = object1_struct;


% -------------------------------- Object 1 -------------------------------

% Object 1 structure
object2_struct = build_obj_struct(CDM_struct.body.object2);

% Reference frame
object2_struct.reference_frame = CDM_struct.body.object2.metadata.contents.REF_FRAME;

% Save in the Conjunction structure
conjunction_struct.object2 = object2_struct;

end


function object_struct_output = build_obj_struct(object_struct_input)

% Quantities which are not essential in the analysis
try
    
    % Maneuverability
    object_struct_output.maneuverable = object_struct_input.metadata.contents.MANEUVERABLE;
    
    % Generic metadata
    object_struct_output.metadata.object_designator = object_struct_input.metadata.contents.INTERNATIONAL_DESIGNATOR;
    object_struct_output.metadata.catalog_name = object_struct_input.metadata.contents.CATALOG_NAME;
    object_struct_output.metadata.object_name = object_struct_input.metadata.contents.OBJECT_NAME;
    object_struct_output.metadata.international_designator = object_struct_input.metadata.contents.INTERNATIONAL_DESIGNATOR;
    object_struct_output.metadata.ephemeris_name = object_struct_input.metadata.contents.EPHEMERIS_NAME;
    object_struct_output.metadata.covariance_method = object_struct_input.metadata.contents.COVARIANCE_METHOD;
catch
end

% Data - radius
area_pc = str2double(object_struct_input.data.additionalParameters.contents.AREA_PC);
object_struct_output.data.radius = sqrt(area_pc/pi)*1e-3; % [km]

% Area-over-mass ratio * Cd
if isfield(object_struct_input.data.additionalParameters.contents,'CD_AREA_OVER_MASS')
    object_struct_output.data.cd_area_over_mass=str2double(object_struct_input.data.additionalParameters.contents.CD_AREA_OVER_MASS);
end

% Area-over-mass ratio * CR
if isfield(object_struct_input.data.additionalParameters.contents,'CR_AREA_OVER_MASS')
    object_struct_output.data.cr_area_over_mass=str2double(object_struct_input.data.additionalParameters.contents.CR_AREA_OVER_MASS);
end

% Data - mean state (km and km/s)
mean_state_object1(1) = str2double(object_struct_input.data.stateVector.contents.X);
mean_state_object1(2) = str2double(object_struct_input.data.stateVector.contents.Y);
mean_state_object1(3) = str2double(object_struct_input.data.stateVector.contents.Z);
mean_state_object1(4) = str2double(object_struct_input.data.stateVector.contents.X_DOT);
mean_state_object1(5) = str2double(object_struct_input.data.stateVector.contents.Y_DOT);
mean_state_object1(6) = str2double(object_struct_input.data.stateVector.contents.Z_DOT);

object_struct_output.data.mean_state = mean_state_object1';


% Data - covariance
covariance_object1(1,1) = str2double(object_struct_input.data.covarianceMatrix.contents.CR_R);

covariance_object1(2,1) = str2double(object_struct_input.data.covarianceMatrix.contents.CT_R);
covariance_object1(2,2) = str2double(object_struct_input.data.covarianceMatrix.contents.CT_T);

covariance_object1(3,1) = str2double(object_struct_input.data.covarianceMatrix.contents.CN_R);
covariance_object1(3,2) = str2double(object_struct_input.data.covarianceMatrix.contents.CN_T);
covariance_object1(3,3) = str2double(object_struct_input.data.covarianceMatrix.contents.CN_N);

covariance_object1(4,1) = str2double(object_struct_input.data.covarianceMatrix.contents.CRDOT_R);
covariance_object1(4,2) = str2double(object_struct_input.data.covarianceMatrix.contents.CRDOT_T);
covariance_object1(4,3) = str2double(object_struct_input.data.covarianceMatrix.contents.CRDOT_N);
covariance_object1(4,4) = str2double(object_struct_input.data.covarianceMatrix.contents.CRDOT_RDOT);

covariance_object1(5,1) = str2double(object_struct_input.data.covarianceMatrix.contents.CTDOT_R);
covariance_object1(5,2) = str2double(object_struct_input.data.covarianceMatrix.contents.CTDOT_T);
covariance_object1(5,3) = str2double(object_struct_input.data.covarianceMatrix.contents.CTDOT_N);
covariance_object1(5,4) = str2double(object_struct_input.data.covarianceMatrix.contents.CTDOT_RDOT);
covariance_object1(5,5) = str2double(object_struct_input.data.covarianceMatrix.contents.CTDOT_TDOT);

covariance_object1(6,1) = str2double(object_struct_input.data.covarianceMatrix.contents.CNDOT_R);
covariance_object1(6,2) = str2double(object_struct_input.data.covarianceMatrix.contents.CNDOT_T);
covariance_object1(6,3) = str2double(object_struct_input.data.covarianceMatrix.contents.CNDOT_N);
covariance_object1(6,4) = str2double(object_struct_input.data.covarianceMatrix.contents.CNDOT_RDOT);
covariance_object1(6,5) = str2double(object_struct_input.data.covarianceMatrix.contents.CNDOT_TDOT);
covariance_object1(6,6) = str2double(object_struct_input.data.covarianceMatrix.contents.CNDOT_NDOT);


covariance_object1 = covariance_object1 + tril(covariance_object1,-1)';

% From m2/s to km2/s -> 1 m2/s = 1*1e-6 km2/s
covariance_object1 = covariance_object1*1e-6;

object_struct_output.data.covariance = covariance_object1;



end


