function opm = opm_fix(opm,m,n)
% Opm_fix writes a structure where the relevant parameters extracted from
% the opm are stored
%
% INPUT: opm = OPM structure from Andrea's function (read_opm)
%        m   = Number of original opm nodes grouped under an optimization
%              node, used to reduce the complexity of the optimization
%        n   = Ending node of the original opm for the simulation
%
% OUTPUT: opm = New OPM structure
%
% Author: Zeno Pavanello, 2023
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

Nfull = length(opm.man_epoch_ignition(1:n));                               % [-] total number of nodes in original opm up to ending point selected by user
Lsc = 6378.173;                                                            % [m] distance scaling factor
Vsc = sqrt(398600.4418/Lsc);                                               % [m/s] velocity scaling factor
Tsc = Lsc/Vsc;                                                             % [s] time scaling factor
Asc = Vsc/Tsc;                                                             % [m/s^2] acceleration scaling factor
dt    = opm.man_duration(1)*m/Tsc;                                         % [-] DeltaT in original opm
epoch_splitted = split(opm.epoch,{'T','-',':'});                          
opm.epoch_spl = epoch_splitted([3,2,1,4,5,6])';                            % [utc] original epoch of the start of the OPM


man_epoch_ignition_spl = split(opm.man_epoch_ignition',{' ','-',':'});   
man_epoch_ignition_spl = man_epoch_ignition_spl(:,[3,2,1,4,5,6]);          % [utc] original epoch of each firing

%% Add nodes in which there is no firing in full time discretization
addNodes = nan(Nfull,1);
et       = nan(Nfull,1);
for i = 1:Nfull-1
    et(i) = utc2et(man_epoch_ignition_spl(i,:))/Tsc;                       % [-] ephemeris time of each fire starting
    if i > 1
        Dt = round(et(i) - et(i-1));
        if Dt > opm.man_duration(i-1)/Tsc                                  % check if fire is continuous wrt last firing
            addNodes(i) = Dt/dt*m;                                         % how many nodes need adding for each original node
        end
    end
end
et(end) = et(end-1) + dt/m;                                                % [-] ephemeris time of the last node
b = round(addNodes(~isnan(addNodes)));                                     % number of nodes to add for each index
a = find(~isnan(addNodes));                                                % indexes in which to add nodes
Nnew = Nfull + sum(b);                                                     % new number of nodes
N    = floor(Nnew/m)+1;                                                    % divide the number of nodes by the required grouping factor m
insert = @(a, x, n)cat(1,  x(1:n), a, x(n+1:end));                         % define insert function
u_rtn = opm.man_dv./opm.man_duration;                                      % [m/s^2] node-wise acceleration of the original OPM
for i = 1:length(b)                                                        % add nodes in which there is no firing
    ind     = a(i) + sum(b(1:i-1)) - 1;
    etAdded = et(ind) + linspace(1,b(i),b(i))'*dt/m;
    et      = insert(etAdded,et,ind);
    uAdded = zeros(3,b(i));
    u_rtn  = cat(2,  u_rtn(:,1:ind), uAdded, u_rtn(:,ind+1:end));
end

%% Create clusters to define new reduced number of nodes
opm.et_start = utc2et(opm.epoch_spl)/Tsc;                                  % [-] ephemeris time of starting of the OPM
uRtn         = nan(3,N);
etNew        = nan(N,1);
for i = 1:N-1                                                              % create clusters of accelerations as mean value inside the interval of length m
    cluster   = m*(i-1)+1:m*i;
    uRtn(:,i) = mean(u_rtn(:,cluster),2)/Asc;
    etNew(i)  = et(cluster(1));                                            % New starting time of each clustered firing
end
etNew(end) = etNew(end-1) + dt;
%constant acceleration between nodes
opm.state    = opm.state./[Lsc*ones(3,1); Vsc*ones(3,1)];                  % [-] Initial state of the primary in cartesian coordinates
states       = nan(6,N);
u            = nan(3,N);
states(:,1)  = opm.state;
for i = 1:N-1                                                              % transform the accelerations into the ECI frame
    u(:,i)  = lvlh2eci(states(1:3,i),states(4:6,i))*lvlh2rtn'*uRtn(:,i);   
    states(:,i:i+1) = propKepOde(states(:,i),u(:,i),[0 dt],1);             % need a propagation to define the transform from RTN to ECI in each node
end
u(:,end)       = zeros(3,1);

%% Define the output
opm.uMax       = max(normOfVec(u));                                        % [m/s^2] Maximum acceleration
opm.u          = u/opm.uMax;                                               % [-] Acceleration in ECI coordinates
opm.uRtn       = uRtn/opm.uMax;                                            % [-] Acceleration in RTN coordinates
opm.nodeStates = states;
opm.finalState = states(:,end);
opm.et         = etNew;
opm.N          = N;                                                        % [-] Total number of nodes for optimization
opm.dt         = dt;                                                       % [-] DeltaT for optimization
opm.canFire    = normOfVec(u) > 0;                                         % Define nodes in which the s/c can fire
opm.coe        = cartesian2kepler(opm.state,1);                            % Initial orbital elements    
opm.T          = 2*pi/opm.coe.n;                                           % [-] Period of the starting orbit
end

