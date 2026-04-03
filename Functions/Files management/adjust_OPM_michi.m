function OPM_struct = adjust_OPM_michi(OPM_struct, sigma_a, sigma_deg )

% DESCRIPTION:
% Adjust the OPM struct for CAM adaptation
% INPUTS:
% - OPM_struct: OPM struct (see read_OPM.m for its initialization)
% - sigma_a: amplitude error in % of maneuver amplitude
%   - size: [1×1 double]
% - sigma_deg: pointing error
%   - size: [1×1 double]
%   - unit: [deg]
% OUTPUT:
% - OPM_struct : OPM struct (see read_OPM.m for its initialization) with
%                additional fields:
%                - man_acc_propagator: array containing the acceleration
%                                      vector in RTN frame
%                  - size: [3×n double]
%                  - unit: [km/s^2]
%                - man_time_propagator: array containing the relevant
%                                       time steps for the maneuver
%                  - size: [1×n double]
%                  - unit: [s]
%                - man_d_time_propagator: array containing the maneuvers
%                                         duration
%                  - size: [3×n double]
%                  - unit: [s]
%                - man_acc_propagator: array containing the covariance
%                                      matrix for simulating process noise
%                  - size: [3×3xn double]
%                  - unit: [km^2/s^4]
%
% Author:  Andrea De Vittori, Politecnico di Milano, 10 April 2022
%          e-mail: andrea.devittori@polimi.it

%  Initialization
size_dt_man = length(OPM_struct.man_epoch_ignition);
time_man_init = zeros(1, size_dt_man);
OPM_struct.man_acc_propagator= [];    % Acceleration array associated with the OPM
OPM_struct.man_time_propagator = [] ; %  Propagation time for the STM method
OPM_struct.man_d_time_propagator = []; % Maneuvering time for the STM method
OPM_struct.man_cov_propagator = [] ;  %  Maneuver covariance for process noise

% Save initial maneuvering time
for i = 1:size_dt_man
    time_man_init(i) = cspice_unitim(cspice_str2et(char(OPM_struct.man_epoch_ignition(i))), 'ET', 'TDT');
end

% Conversion from ET to UTC
time_man_init = time_man_init - cspice_unitim(cspice_str2et(char(OPM_struct.epoch)), 'ET', 'TDT');

% check for null maneuver at initial time:
if time_man_init(1)>0
    OPM_struct.man_time_propagator(end+1) = 0;
    OPM_struct.man_d_time_propagator(end+1) = time_man_init(1);
    OPM_struct.man_acc_propagator(:, end+1) =  zeros(3, 1);%# RTN frame
    % The acceleration Covariance is null for non thrusting arcs
    OPM_struct.man_cov_propagator(:,:,end) = zeros(3, 3); %# RTN frame
elseif time_man_init(1)<0
    manindex=find(time_man_init>=0);
    if time_man_init(manindex(1))==0
        OPM_struct.man_dv(:,1:manindex(1)-1)=[];
        OPM_struct.man_epoch_ignition(1:manindex(1)-1)=[];
        OPM_struct.man_duration(1:manindex(1)-1)=[];
        size_dt_man = length(OPM_struct.man_epoch_ignition);
        time_man_init = time_man_init(manindex(1):end);
    elseif time_man_init(manindex(1))>0
        dt=man_duration(mainindex(1)-1);
        dt2=time_man_init(manindex(1));
        maneuver_tbd=OPM_struct.man_dv(:,manindex(1)-1)*dt2/dt;
        OPM_struct.man_dv(:,1:manindex(1)-1)=[];
        OPM_struct.man_epoch_ignition(1:manindex(1)-1)=[];
        OPM_struct.man_duration(1:manindex(1)-1)=[];

        OPM_struct.man_dv=[maneuver_tbd,OPM_struct.man_dv];
        OPM_struct.man_epoch_ignition=[OPM_struct.epoch, OPM_struct.man_epoch_ignition];
        OPM_struct.man_duration=[dt2,OPM_struct.man_duration];

        size_dt_man = length(OPM_struct.man_epoch_ignition);
        time_man_init = [0, time_man_init(manindex(1):end)];
    end

end

% Adjust OPM struct for CAM adaptation
for i= 1:size_dt_man
    % Saving info when the thruster is on
    OPM_struct.man_time_propagator(end+1) = time_man_init(i);
    OPM_struct.man_d_time_propagator(end+1) =  OPM_struct.man_duration(i);
    OPM_struct.man_acc_propagator(:, end+1) = OPM_struct.man_dv(:, i)/OPM_struct.man_duration(i);
    intermediate_time = time_man_init(i) + OPM_struct.man_duration(i);
    tau = reshape(OPM_struct.man_acc_propagator(:,end)/norm(OPM_struct.man_acc_propagator(:,end)), [3,1]);

    % Build the Acceleration pointing and magnitude covariance matrix
    a_T = norm(OPM_struct.man_acc_propagator(:, end));
    [~, alpha, beta] = cspice_reclat(OPM_struct.man_acc_propagator(:, end));
    S = [-cos(beta)*sin(alpha), -cos(alpha)*sin(beta); ...
        cos(alpha)*cos(beta), -sin(alpha)*sin(beta);...
        0, cos(beta)];
    EPS = eye((2))*((sigma_deg*pi/180.*1./3.)^2);
    R = (a_T*sigma_a/100.*1./3.)^2*dot(tau, tau) + a_T^2* S*EPS*S';
    if isempty (OPM_struct.man_cov_propagator)
        OPM_struct.man_cov_propagator(:,:, end) = R;% RTN frame;
    else
        OPM_struct.man_cov_propagator(:,:, end+1) = R;% RTN frame;
    end

    % Saving info when the thruster is off
    if i <= size_dt_man - 1
        if intermediate_time < time_man_init(i+1)
            OPM_struct.man_time_propagator(end+1) = intermediate_time;
            OPM_struct.man_d_time_propagator(end+1) = time_man_init(i+1) - intermediate_time;
            OPM_struct.man_acc_propagator(:, end+1) =  zeros(3, 1);%# RTN frame
            % The acceleration Covariance is null for non thrusting arcs
            OPM_struct.man_cov_propagator(:,:,end+1) = zeros(3, 3); %# RTN frame
        end
    end
end
end

