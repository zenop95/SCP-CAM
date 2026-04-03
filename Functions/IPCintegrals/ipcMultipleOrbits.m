function [err, maxIpcConvolutes, ipc] = ipcMultipleOrbits(sample,N,varargin)
%ipcMultipleOrbits Runs the caseStudy script for different values of the
%initial conditions.

%% Initialize variables

%gravitational constant
mu = 3.986004418e14; %[m^3/s^2]

%Hard Body Radius
HBR = 100; %[m] 

% primary orbit
a      = sample.a;      %[m] semi-major axis
e      = sample.e;      %[-] eccentricity
theta0 = sample.TA;     %[rad] initial true anomaly
n      = sqrt(mu/a^3);   %[rad/s] mean motion
T      = 2*pi/n;         %[s] orbital period

% Relative orbit
C0        = sample.C0;
expValue0 = sample.mu0';

% time grid
t0     = 0;       %[s] initial time
tf     = 2*T;     %[s] final time
t      = linspace(t0,tf,N);  

STMGrid         = nan(6,6,N);
C               = nan(6,6,N);
posC            = nan(3,3,N);
expValue        = nan(6,N);
posExpValue     = nan(3,N);
normPosExpValue = nan(1,N);
cubeIpc         = nan(1,N);
sphereIpc       = nan(1,N);
numIpc          = nan(1,N);
constIpc        = nan(1,N);
xyIpc           = nan(1,N);
xzIpc           = nan(1,N);
yzIpc           = nan(1,N);
maxIpc          = nan(1,N);

%% Analytical/numerical methods
for i =1:N
% compute YA State transition matrix, mean and covariance of relative position
% at each point    
    if e > 0.001
        STM = YAStateTransition(a,e,theta0,t(i),t0,mu);
    else
        STM = CWStateTransition(a,t(i),t0,mu);
    end
    STMGrid(:,:,i)     = STM;
    C(:,:,i)           = STM*C0*STM';
    posC(:,:,i)        = C(1:3,1:3,i);
    expValue(:,i)      = STM*expValue0;
    posExpValue(:,i)   = expValue(1:3,i);
    normPosExpValue(i) = norm(posExpValue(:,i));

    % Compute all the required IPC integrals
    [maxIpc(i),cubeIpc(i),sphereIpc(i),constIpc(i),numIpc(i),xyIpc(i),...
        xzIpc(i),yzIpc(i)] = ipcWrapper(posExpValue(:,i),posC(:,:,i), ...
        HBR,'max','cuboid','spheroid','constant','precise','xy','xz','yz');
end
ipc = struct( ...
            'max', maxIpc, ...
            'cuboid', cubeIpc, ...
            'spheroid', sphereIpc, ...
            'constant', constIpc, ...
            'precise', numIpc, ...
            'xy', xyIpc, ...
            'xz', xzIpc, ...
            'yz', yzIpc ...
            );

% Compute absolute and relative error for each integral performed
[maxErr,cubeErr,sphereErr,constErr,xyErr,xzErr,yzErr] = ...
    ipcErrors(numIpc,maxIpc,cubeIpc,sphereIpc,constIpc,xyIpc,xzIpc,yzIpc);

err = struct( ...
            'maxErr', maxErr, ...
            'cubeErr', cubeErr, ...
            'sphereErr', sphereErr, ...
            'constErr', constErr, ...
            'xyErr', xyErr, ...
            'xzErr', xzErr, ...
            'yzErr', yzErr ...
            );
a = maxIpc > numIpc;
maxIpcConvolutes.pointwise = a;
b = 0;
for i = 1:N
    b = b + double(a(i) == 0);
end
maxIpcConvolutes.total = b/N; %percentage of time for which ICP_max<ICP 
%% plots
if nargin > 2 && strcmpi(varargin{1},'plots')

    figure(1)
    plot(t/T,normPosExpValue)
    xlabel('t/T_{orb} [-]')
    ylabel('\mu_r [m]')
    grid on
    
    figure(2)
    semilogy(t/T,maxIpc)
    hold on
    semilogy(t/T,cubeIpc)
    semilogy(t/T,sphereIpc)
    semilogy(t/T,constIpc)
    semilogy(t/T,numIpc,'k')
    % semilogy(mcT/T,mcIpc,'k*')
    xlabel('t/T_{orb} [-]')
    ylabel('IPC [-]')
    legend('IPC_{max}','Cuboid IPC','Spheroid IPC','Constant IPC','Precise IPC')
    ylim([5e-8, 1])
    grid on
    hold off
    
    figure(3)
    semilogy(t/T,xyIpc)
    hold on
    semilogy(t/T,xzIpc)
    semilogy(t/T,yzIpc)
    semilogy(t/T,numIpc,'k')
    xlabel('t/T_{orb} [-]')
    ylabel('IPC [-]')
    legend('XY projected IPC','XZ projected IPC','YZ projected IPC','Precise IPC')
    ylim([5e-8, 1])
    grid on
    hold off
end