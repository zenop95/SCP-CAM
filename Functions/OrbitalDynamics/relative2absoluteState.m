function orbit = relative2absoluteState(chaser, target, varargin)
% The function provides cartesian vector r (position) and v (velocity) of
% the Chaser by relative position wrt the Target.
%
%   Usage:
%       orbit = orbitFromRelativePosition(chaser, target,...
%           'Position',                 [0 0 0], ...
%         	'Velocity',                 [0 0 0], ...
%           'PositionFrame',            'LVLH', ...
%           'Attitude',                 [0 pi 0] || [1 0 0 0] || eye(3), ...
%           'AimedAxis',                '+Z', ...
%           'AttitudeFrame',            'TAPF', ...
%           'ChaserAttitudeFrame',      'CROF', ...
%           'EulAngRotationSequence',   'zyx', ...
%           'AttitudeMatrixType',       'dcm', ...
%           'MaxOOPAmplitude',          20);
%
%   Inputs:
%       - target: a structure defining the target geometric properties and
%           orbit.
%
%   Parameters:
%       - Position: position, in [m], of the chaser wrt the target,
%           expressed in the later specified frame; 3x1 vector.
%           Default: [0 0 0].
%       - PositionFrame: the frame with respect to which the chaser
%           position shall be computed.
%           Accepted position frames: 'LVLH' | 'TROF' | 'TGFF' | 'TAPF'
%           Default: 'LVLH'.
%       - Attitude: attitude of the chaser; destination and origin frames
%           are specified by a latter parameter.
%           Accepted formats are: Euler angles (3x1), Quaternions (4x1) or
%           DCM / Rotation matrix (3x3); the type of the latter can be
%           specified by a dedicated parameter.
%       - AimedAxis: alternative to Attitude; specifies which axis of the
%           ChaserAttitudeFrame to aim at the center of AttitudeFrame. In
%           this configuration, the sole supported AttitudeFrame are 'LVLH'
%           and 'TROF' - which are concident.
%           Accepted aimed axis: '±X' | '±Y' | '±Z'
%           Default: '' (i.e., aimed axis not considered).
%       - AttitudeFrame: the frame in which the chaser's attitude shall be
%           resolved.
%           Accepted reference attitude frames: 'ECI' | 'LVLH' | 'TROF' |
%           'TGFF' | 'TAPF'.
%           Note:   'TAPF' only supported in AimedAxis mode
%                   'ECI' only in Attitude mode.
%           Default: 'LVLH'.
%       - ChaserAttitudeFrame: the frame of the chaser required to be
%           rotated wrt the provided AttitudeFrame. E.g., by specifying
%           'CSFF' as ChaserAttitudeFrame, 'LVLH' as AttitudeFrame and
%           [0 0 0] as Euler sequence, the 'CSFF' will be aligned with the
%           orbital frame LVLH at the start of the simulation.
%           The specified frame will be neglected if AimedAxis is used.
%           Accepted chaser attitude frames: 'CROF' | 'CSFF' | 'CGFF'.
%           Default: 'CROF'.
%       - EulAngRotationSequence: specifies the rotation sequence to be
%           used if a 3x1 vector is provided as Attitude parameter.
%           The parameter is neglected if a Euler sequence is not used.
%           Default: 'zyx'.
%       - AttitudeMatrixType: specifies whether the provided 3x3 matrix is
%           a Direction Cosine or Rotation matrix.
%           Accepted matrix: 'DCM', 'ROTMAT'.
%           Default: 'DCM'.
%
%   Output:
%       - orbit: a structure containing the descriptive parameters of the
%           chaser's orbit in terms of orbital radius [km] and
%           velocity [km/s]. Both are expressed in ECI.

%% Inputs parsing
% Parse and validate the inputs
parser = inputParser;

% Required objects
addRequired(parser, 'chaser');
addRequired(parser, 'target');

% Relative/Absolute positioning (3x1);
addParameter(parser, 'Position', zeros(3, 1), ...
    @(x) validateattributes(x, {'numeric'}, {'real', 'finite', 'nonnan', 'size', [3, 1]}));
% Positioning frame name (string)
addParameter(parser, 'PositionFrame', 'LVLH', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
% Relative/Absolute velocity (3x1);
addParameter(parser, 'VelType', 'lvlhSyncd', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(parser, 'Velocity', zeros(3, 1), ...
    @(x) validateattributes(x, {'numeric'}, {'real', 'finite', 'nonnan', 'size', [3, 1]}));
% Relative/Absolute attitude (variable)
% Defaults to the chaser body (CROF) being parallel to LVLH
addParameter(parser, 'Attitude', [0 0 0], ...
    @(x) validateattributes(x, {'numeric'}, {'real', 'finite', 'nonnan'}));
% Axis aiming (pointing attitude)
addParameter(parser, 'AimedAxis', '', ...
    @(x) validateattributes(x, {'char', 'string'}, {}));
% Attitude resolution frame
addParameter(parser, 'AttitudeFrame', 'LVLH', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
% Chaser attitude definition frame
addParameter(parser, 'ChaserAttitudeFrame', 'CROF', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Specify the rotation sequence if a 3x1 or 1x3 vector or 3x3 matrix is
% provided as chaser attitude
addParameter(parser, 'EulAngRotationSequence', 'zyx', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
% Specify the type of the rotation matrix if a 3x3 variable is provided as
% chaser attitude (DCM | ROTMAT)
addParameter(parser, 'AttitudeMatrixType', 'dcm', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Additional parameter for the out-of-plane motion; if set <0 the provided
% starting condition will be used as maximum Y displacement.
addParameter(parser, 'MaxOOPAmplitude', -1, ...
    @(x) validateattributes(x, {'numeric'}, {'real', 'finite', 'nonnan', 'scalar'}));

% Relative/Absolute body rates (3x1); defaults to [0 0 0] in CROF w.r.t.
% TROF
addParameter(parser, 'BodyRates', zeros(3, 1), ...
    @(x) validateattributes(x, {'numeric'}, {'real', 'finite', 'nonnan', 'size', [3, 1]}));
%Specify the reference frame with respect to which the angular rates of the
%chaser are expressed
addParameter(parser, 'BodyRatesFrame', 'TROF', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

parse(parser, chaser, target, varargin{:});

% Accepted inputs
validPosFrames = {'lvlh', 'trof', 'tapf', 'tgff'};
validAttFrames = {'eci', 'lvlh', 'trof', 'tgff', 'tapf'};
validChaserAttFrames = {'csff', 'cgff', 'crof'};
validRatesFrames = {'trof','lvlh','eci'};
validAimedAxis = {'', '+x', '+y', '+z', '-x', '-y', '-z'};
validVelType = {'lvlhSyncd', 'trofSyncd', 'custom'};

% Position frame validation
posFrame = lower(parser.Results.PositionFrame);
if ~any(strcmp(validPosFrames, posFrame))
    error('Invalid position frame.');
end
% Attitude frame validation
attDestFrame = lower(parser.Results.AttitudeFrame);
if ~any(strcmp(validAttFrames, attDestFrame))
    error('Invalid reference attitude frame.');
end
% Chaser attitude frame validation
chAttFrame = lower(parser.Results.ChaserAttitudeFrame);
if ~any(strcmp(validChaserAttFrames, chAttFrame))
    error('Invalid chaser attitude frame.');
end
% Rates frame validation
ratesFrame = lower(parser.Results.BodyRatesFrame);
if ~any(strcmp(validRatesFrames, ratesFrame)) 
    error('Invalid body rates frame.');
end
% Aimed axis validation
aimedAxis = lower(parser.Results.AimedAxis);
if ~any(strcmp(validAimedAxis, aimedAxis))
    if ~any(strcmp({'+','-'}, aimedAxis(1)))
        error('You must specify the sign (±) before the axis. + shall not be neglected');
    else
        error('Invalid aimed axis definition..');
    end
end

velType = parser.Results.VelType;
if ~any(strcmp(validVelType, velType))
    error('Invalid velocity type.');
end

% Relative distance computed between specified frame and CROF
% (Chaser's centers of mass)
% Target's position and velocity in ECI
tr = target.orbit.eci2trof.translation;
tv = target.orbit.eci2trof.velocity;

% DCM matrix rotating LVLH into ECI and angular rate of ECI wrt LVLH in
% LVLH
[DCM_IO, w_IO_O] = OrbitalDynamics.eci2lvlh(tr, tv);
nn = norm(w_IO_O);

DCM_TI = target.orbit.eci2trof.dcm;
DCM_TO = DCM_TI * DCM_IO;

% Reference systems stransform
tgff2trof = target.transforms.tgff2trof;
tgff2tapf = target.transforms.tgff2tapf;
crof2csff = chaser.transforms.crof2csff;
cgff2crof = chaser.transforms.cgff2crof;
cgcf2crof = chaser.transforms.cgcf2crof;


%% Attitude
% Nomenclature:
% - Chaser: specified "from" frame - chAttFrame
% - CROF: chaser's body frame
% - LVLH: orbital reference frame
% - TROF: target's body frame
% - Destination: specified "to" frame - attDestFrame
% - ECI: inertial frame

% If no aligned axis is provided, use absolute attitude definition
% Raw unprocessed input
chaserAttitudeRaw = parser.Results.Attitude;
attitudeSize = size(chaserAttitudeRaw);

% Converts the input type to a DCM
if isequal(attitudeSize, [4, 1]) || isequal(attitudeSize, [1, 4])
    % Quaternion
    dcmChaserToDest = quat2dcm(chaserAttitudeRaw);
elseif isequal(attitudeSize, [3, 1]) || isequal(attitudeSize, [1, 3])
    % Eulers angles; rotation sequence provided as input parameter
    seq = parser.Results.EulAngRotationSequence;
    dcmChaserToDest = eulang2dcm(reshape(chaserAttitudeRaw,3,[]), seq);
elseif isequal(attitudeSize, [3, 3])
    % Check whether it's a matrix or a DCM
    if strcmpi(parser.Results.ChaserAttFrame, 'rotmat')
        dcmChaserToDest = chaserAttitudeRaw;    % Rotation matrix
    elseif strcmpi(parser.Results.ChaserAttFrame, 'dcm')
        dcmChaserToDest = chaserAttitudeRaw';   % DCM
    else
        error('Unknown matrix type.');
    end
else
    error('Unknown attitude input size.');
end

% Destination attitude frame; rotating DEST into LVLH
switch attDestFrame
    case 'eci'
        dcmDestToLvlh = DCM_IO';
    case {'trof', 'tgff'}
        % TGFF, TAPF and TROF have aligned axis -> same case
        % target.orbit.attitude.rotmat' is
        dcmDestToLvlh = DCM_TI * DCM_IO;
    case 'lvlh'
        dcmDestToLvlh = eye(3);
    otherwise
        error('%s frame not supported for attitude mode Attitude', upper(attDestFrame));
end

% Chaser attitude definition frame; rotating the provided frame into CROF
switch chAttFrame
    case {'cgff', 'crof'}
        % CGFF and CROF have aligned axis -> same case
        dcmChaserToCrof = eye(3);
    case 'csff'
        dcmChaserToCrof = chaser.body.cgff2crof.dcm*chaser.body.crof2csff.dcm;
        % DCM_bodytoCSFF'
    otherwise
        error('Unknown Chaser attitude frame.');
end

% Rotation of the Chaser is composed of the rotation from the required
% Chaser definition frame to the CROF one, then from CROF to LVLH, finally
% from LVLH to ECI. Attitude is defined in ECI.
%     dcmChaserEciToBody = dcmChaserToCrof * dcmChaserToDest' * ...
%         dcmDestToLvlh' * DCM_IO';
DCM_CI = dcmChaserToCrof' * dcmChaserToDest * ...
dcmDestToLvlh * DCM_IO';
% Providing the attitude rotation matrix
orbit.eci2crof.dcm = DCM_CI;

DCM_CT = DCM_CI*DCM_TI';
%% Relative position (LVLH)
% Define rotation and position offset
% - posOffset: vector from TROF to ref. frame, in TROF axis
switch posFrame
    case 'trof'
        DCM_TP    = eye(3);
        posOffset = zeros(3, 1);
    case 'tapf'
        DCM_TP    = tgff2trof.dcm*tgff2tapf.dcm';
        posOffset = tgff2trof.dcm*(tgff2tapf.translation - ...
                                    tgff2trof.translation);
    case 'tgff'
        % Aligned with TROF
        DCM_TP    = tgff2trof.dcm;
        posOffset = -tgff2trof.dcm*tgff2trof.translation;
    case 'lvlh'
        DCM_TP    = DCM_TO;
        posOffset = -target.orbit.lvlh2trof.dcm*...
                     target.orbit.lvlh2trof.translation;
    otherwise
        error('Unexpected positioning frame.');
end
crofOffset = (DCM_CT'*cgcf2crof.dcm)*cgcf2crof.translation;
% Relative position of CROF wrt TROF in TROF
drTrof = DCM_TP * parser.Results.Position + posOffset + crofOffset;

% Relative position of CGCF wrt TROF in LVLH frame
drLvlh = DCM_TO' * drTrof; %[m]

% Absolute position of chaser in ECI frame
orbit.eci2crof.translation = tr + DCM_TI' * drTrof;

%% Angular velocity
%rateFrame is always expressed in CROF wrt the selected frame
bodyRates = parser.Results.BodyRates;
DCM_CO    = DCM_CI*DCM_IO;
% angular rates of ECI wrt LVLH in CROF
w_IO_C    = DCM_CO*w_IO_O;
switch ratesFrame
    case 'trof'
       orbit.eci2crof.bodyRates = bodyRates + ...
           DCM_CT*target.orbit.eci2trof.bodyRates;
       orbit.lvlh2crof.bodyRates = orbit.eci2crof.bodyRates + w_IO_C;
    case 'lvlh'
        % Angular rates of Body wrt LVLH (expressed in body frame) [rad/s]
        orbit.lvlh2crof.bodyRates = bodyRates;
        % Angular rates of Body wrt to ECI (expressed in body frame) [rad/s]
        orbit.eci2crof.bodyRates  = (bodyRates - w_IO_C);
    case 'eci'
        % Angular rates of Body wrt to ECI (expressed in body frame) [rad/s]
        orbit.eci2crof.bodyRates = bodyRates;
        % Angular rates of Body wrt LVLH (expressed in body frame) [rad/s]
        orbit.lvlh2crof.bodyRates = bodyRates + w_IO_C;
    otherwise
        error('%s frame not supported for body rates definition', upper(ratesFrame));
end

%% Absolute velocity of chaser in ECI frame
switch velType
    case 'lvlhSyncd'
        orbit.lvlh2crof.velocity = zeros(3,1);   %[m/s]
    case 'trofSyncd'
        %angular velocity of TROF wrt LVLH in LVLH
        w_TO_O = target.orbit.lvlh2trof.dcm'*...
                        target.orbit.lvlh2trof.bodyRates;
        %relative velocity of CROF wrt LVLH in LVLH
        orbit.lvlh2crof.velocity = cross(w_TO_O, drLvlh); %[m/s]
    case 'custom'
        orbit.lvlh2crof.velocity = parser.Results.Velocity; %[m/s]
    otherwise
        error('Invalid velocity type definition')
end
velLvlh = orbit.lvlh2crof.velocity;
orbit.lvlh2crof.translation = drLvlh + target.orbit.lvlh2trof.translation; %[m]
orbit.lvlh2crof.dcm = DCM_CO;
orbit.eci2crof.velocity = tv + DCM_IO*velLvlh - DCM_IO*cross(w_IO_O, drLvlh);
end