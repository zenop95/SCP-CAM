function [const,max,cube] = ipcDaComparison(sample,t0,t,HBR)
%ipcDaComparison Runs the caseStudy script for different values of the
%initial conditions.

%% Initialize variables

%gravitational constant
mu = 3.986004418e14; %[m^3/s^2]

% primary orbit
a      = sample.a;      %[m] semi-major axis
e      = sample.e;      %[-] eccentricity
theta0 = sample.TA;     %[rad] initial true anomaly

% Relative orbit
C0          = sample.C0;
expValue0   = sample.mu0;
N           = length(t);
C           = nan(6,6,N);
posC        = nan(3,3,N);
expValue    = nan(6,N);
posExpValue = nan(3,N);
const   = nan(N,1);
max     = nan(N,1);
cube     = nan(N,1);

%% Analytical/numerical methods
for i =1:N
% compute YA State transition matrix, mean and covariance of relative position
% at each point    
    if e > 0.001
        STM = YAStateTransition(a,e,theta0,t(i),t0,mu);
    else
        STM = CWStateTransition(a,t(i),t0,mu);
    end
    C(:,:,i)           = STM*C0*STM';
    posC(:,:,i)        = C(1:3,1:3,i);
    expValue(:,i)      = STM*expValue0;
    posExpValue(:,i)   = expValue(1:3,i);
    % Compute all the required IPC integrals
    [const(i,:),max(i,:),cube(i,:)]  = ...
        ipcWrapper(posExpValue(:,i),posC(:,:,i),HBR,'constant','max','cuboid');
end