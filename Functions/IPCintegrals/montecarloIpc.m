function Ipc = montecarloIpc(mu0,P0,R,a,e,theta0,t,t0,M,varargin)
% montecarloIpc Performes a set of montecarlo simulations to estimate the
% Instantaneous probability of collision for a given relative position
% represented by a normally distributed random variable. The highest the
% number of samples, the more precise the accuracy of the estimation.

% INPUT: mu0 [m]      = expected value of the relative position random variable
%        P0 [m.^2]    = covariance Matrix of the relative position random variable
%        R [m]        = hard Body Radius (HBR) of the collision
%        a [m]        = semimajor axis of the orbit of the spacecraft
%        e [-]        = eccentricity of the orbit of the spacecraft
%        theta0 [rad] = offset true anomaly
%        t [s]        = end time of the evolution of the dynamics represented by
%                       the state transition matrix
%        t0 [s]       = start time of the evolution of the dynamics represented by
%                       the state transition matrix
%        M [-]        = number of Montecarlo samples

% OUTPUT: 
%        Ipc = [-] Computed Instantaneous Collision Probability

% Author: Zeno Pavanello 2021
%--------------------------------------------------------------------------

if nargin == 10
    mu = varargin{1}; 
else 
    mu = 398600.4418; 
end

N            = length(t);
state        = nan(6,M);
normDist     = nan(N,M);
mcSamples    = mvnrnd(mu0,P0,M); % montecarlo samples
nContact     = zeros(N,1);
STMGrid      = nan(6,6,N);

%% compute YK State transition matrix, mean and covariance of relative position
% at each point
for i = 1:N
    if e > 0.001
        STM = YAStateTransition(a,e,theta0,t(i),t0,mu);
    else
        STM = CWStateTransition(a,t(i),t0,mu);
    end
    STMGrid(:,:,i) = STM;
end

%% Montecarlo Simulations
for j = 1:M
    for i = 1:N
        state(:,i)    = STMGrid(:,:,i)*mcSamples(j,:)';
        normDist(i,j) = norm(state(1:3,i)); %samples in columns, time in rows
        nContact(i)   = nContact(i) + (normDist(i,j) <= R);
    end
end
Ipc = nContact/M;
end

