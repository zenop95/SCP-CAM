function OPM_struct = read_OPM(OPM_filename)
% Create a structure containing all the OPM file information
% -------------------------------------------------------------------------
% INPUTS:
% - OPM_filename : OPM path file
% -------------------------------------------------------------------------
% OUTPUT:
% - OPM_struct : structure whose fields are:
%               - filename: OPM path file
%               - state: state at initial maneuvering time t0
%                 - size: [6×1 double]
%                 - unit: [km, km/s]
%               - covariance: covariance matrix at initial maneuvering time
%                             t0
%                  - size: [6×6 double]
%                  - unit: [km^2, km^2/s; km^2/s, km^2/s^2]
%               - man_dv: delta v array associated to each maneuver branch
%                 - size:[3×n double]
%                 - unit: [km/s]
%               - man_epoch_ignition: starting date array for each maneuver
%                                     branch
%                  - size:[1×n str]
%               - man_duration: array containing the maneuver duration
%                 - size:[1×n double]
%                 - unit: [s]
%               - epoch: initial maneuver date
%                 - size:[1×1 str]
%                 - unit: [s]
%               - ccsds_opm_vers: OPM version
%                 - size:[1×1 str]
%               - creation date: OPM creation date
%                 - size:[1×1 str]
%               - originator: Company who generated the OPM file
%                 - size:[1×1 str]
%               - object_name: space object name
%                 - size:[1×1 str]
%               - object_id: space object id
%                 - size:[1×1 str]
%               - center_name: reference frame origin
%                 - size:[1×1 str]
%               - ref_frame: reference frame
%                 - size:[1×1 str]
%               - time_system: time reference frame
%                 - size:[1×1 str]
%               - mass: space object mass at t0
%                 - size:[1×1 str]
%                 - unit: [kg]
%               - solar_rad_area: area subjected to SRP
%                 - size:[1×1 str]
%                 - unit: [m^2]
%               - solar_rad_coeff: solar radiation coeffiecient
%                 - size:[1×1 str]
%               - area: sapce object area
%                 - size:[1×1 double]
%                 - unit: [m^2]
%               - Cd: space object drag coefficient
%                 - size:[1×1 double]
%               - cov_ref_frame: covariance reference frame
%                 - size:[1×1 str]
%               - man_delta_mass: fuel mass cost
%                 - size:[1×1 str]
%               - man_ref_frame: maneuver reference frame
%                 - size:[1×1 str]
%               - user_defined_ccd_cd:
%                 - size:[1×1 str]
%               - b_star: ballistic coefficient:
%                 - size:[1×1 double]
% -------------------------------------------------------------------------
%
% Author: Andrea De Vittori, Politecnico di Milano, 10 July 2021
%         e-mail: andrea.devittori@polimi.it


% Open the file
fid = fopen(OPM_filename, 'r');

% Display an error message if the file is not found
if fid==-1
    ME = MException('PoliMiSST:OPMFileError:fileNotFound', ...
        'The OPM file was not found');
    throw(ME)
end

% Key words to look up in the OPM file
keys = ["X", "Y", "Z", "X_DOT", "Y_DOT", "Z_DOT", "CX_X", "CY_X",...
    "CZ_X", "CX_DOT_X","CY_DOT_X", "CZ_DOT_X", "CY_Y", "CZ_Y",...
    "CX_DOT_Y",  "CY_DOT_Y", "CZ_DOT_Y", "CZ_Z", "CX_DOT_Z",...
    "CY_DOT_Z", "CZ_DOT_Z", "CX_DOT_X_DOT", "CY_DOT_X_DOT",...
    "CZ_DOT_X_DOT", "CY_DOT_Y_DOT", "CZ_DOT_Y_DOT", "CZ_DOT_Z_DOT",...
    "MAN_EPOCH_IGNITION", "MAN_DURATION", "MAN_DV_1", "MAN_DV_2", "MAN_DV_3", "EPOCH", "MASS", "DRAG_AREA", "DRAG_COEFF"];

% Saving the position where the covariance key starts
pos_init_cov = find(keys == "CX_X");
pos_time_ignition = find(keys =="MAN_EPOCH_IGNITION");
pos_dv = find(keys =="MAN_DV_1");
pos_epoch = find(keys =="EPOCH");
content = fileread(OPM_filename) ;
occurrences = numel(regexp(content, 'MAN_DV_1'));

% mean_state vector and covariance matrix/el initialization
OPM_struct = struct();
OPM_struct.filename = OPM_filename; % OPM path file
OPM_struct.state = zeros(6, 1); % State at t0 (initial maneuvering point)
covariance_el = zeros(21, 1); % Covariance elements at t0 (it's symmetric and only 21 elements are taken)
covariance_matrix = tril(ones(6));  % Covariance matrix at t0
OPM_struct.covariance = triu(ones(6));% Covariance matrix at t0
OPM_struct.man_dv = zeros(3, occurrences) ; % Delta V for each maneuver arc in RTN
OPM_struct.man_epoch_ignition = strings([1,occurrences]); % Starting maneuvering time for each thrust arc
OPM_struct.man_duration = zeros(1,occurrences); % Thrust arc durations
OPM_struct.epoch = strings([1,1]); % Epoch corresponding to t0

% Line initialization for showing the file line where an error occurs
n_line = 0;
index = 0;

% Reading the file
while ~feof(fid)
    % Line conversion to a string
    line = string(strip(fgetl(fid)));
    n_line = n_line + 1;
    % Skip empty lines
    if line == ""
        continue
    end
    line = strip(split(line, "="));

    % Check the correct line length
    if size(line, 1) == 2
        pos_key = find(keys == line(1));
        if isempty(pos_key) % not mean state nor covariance
            % Save the field in the struct
            OPM_struct.(lower(line(1))) = line(2);
        elseif pos_key < pos_init_cov
            OPM_struct.state(pos_key) = str2double(line(2));
        elseif  pos_key < pos_time_ignition
            % Temporarily store covariance matrix as 1D array
            covariance_el(pos_key - pos_init_cov +1) = str2double(line(2));
        elseif strcmp(line(1), 'MAN_EPOCH_IGNITION')
            index = index+1;
            OPM_struct.(lower(line(1)))(index) = erase(line(2), 'Z');
            OPM_struct.(lower(line(1)))(index) = replace(OPM_struct.(lower(line(1)))(index), 'T', ' ');
        elseif strcmp(line(1),'MAN_DURATION')
            % save maneuver duration
            OPM_struct.(lower(line(1)))(index) = double(line(2));
        elseif pos_key >= pos_dv && pos_key < pos_epoch
            % save maneuver delta v
            OPM_struct.man_dv(floor(str2double(erase(line(1),'MAN_DV_'))), index) = str2double(line(2));
        elseif pos_key == pos_epoch
            % correct data in utc format
            OPM_struct.(lower(line(1)))(1) = erase(line(2), 'Z');
        elseif strcmp(line(1), 'MASS')
            % save mass
            OPM_struct.mass =  double(line(2));
        elseif strcmp(line(1), 'DRAG_AREA')
            % save area
            OPM_struct.area =  double(line(2));
        elseif strcmp(line(1), 'DRAG_COEFF')
            % save drag coefficient
            OPM_struct.Cd   =  double(line(2));
        end
    elseif startsWith(line, "COMMENT")
        % Collect comments
        try
            k = size(OPM_struct.comment, 2);
            OPM_struct.comment(k+1) = strip(extractAfter(line, "COMMENT"));
        catch
            OPM_struct.comment(1) = strip(extractAfter(line, "COMMENT"));
        end
    else
        % Raise an error if the OPM is not in compliance with the standard
        ME = MException("PoliMiSST:OPMFileError:ErrorReadData", ...
            "The OPM file is not consistent with the standard at line: %d", ...
            n_line);
        throw(ME)
    end
end

% save bstar
OPM_struct.bstar = OPM_struct.area/OPM_struct.mass*OPM_struct.Cd;

% Save the mean_state and covariance matrix in the struct
if any(covariance_el(:))
    covariance_matrix(covariance_matrix ~= 0) = covariance_el;
    OPM_struct.covariance = covariance_matrix + tril(covariance_matrix, -1).';
end

% Close the file
fclose(fid);
end
