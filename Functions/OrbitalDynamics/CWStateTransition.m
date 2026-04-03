function CWMatrix = CWStateTransition(a,t,t0,varargin)
%CWSTATETRANSITION outputs the Clohessy-Wiltshire state transition matrix
%that can be used to compute the linearized dynamics of the relative motion
%between to spacecraft in LVLH frame. 
% INPUT: a [km]     = mean motion of the orbit of the spacecraft
%        t [s]      = end time of the evolution of the dynamics represented by
%                     the YA state transition matrix
%        t0 [s]     = start time of the evolution of the dynamics represented by
%                     the YA state transition matrix
%        mu         = standard gravitational constant of the central body
%        refFrame   = (char) Reference frame in which the input state
%                       vector is expressed. Can either be 'Hill' or 'Lvlh'.

% OUTPUT: CWMatrix = Clohessy-Wiltshire state transition matrix
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

if nargin > 3 
    mu = varargin{1}; 
else 
    mu = 398600.4418; 
end

Dt  = t-t0;
n   = sqrt(mu/a^3);
c   = cos(n*Dt);
s   = sin(n*Dt);

CWMatrix = [4-3*c 0 0 s/n 2*(1-c)/n 0;
           6*(s-n*Dt) 1 0 2*(c-1)/n 4*s/n-3*Dt 0;
           0 0 c 0 0 s/n;
           3*n*s 0 0 c 2*s 0;
           6*(c-1)*n 0 0 -2*s 4*c-3 0;
           0 0 -n*s 0 0 c];

if nargin > 4 
    refFrame = lower(varargin{2});
    dcm       = lvlh2hill()';
    Hill2Lvlh = [dcm zeros(3,3); zeros(3,3) dcm];
    switch refFrame
    case 'rtn'
        return;
    case 'lvlh'
        CWMatrix = Hill2Lvlh*CWMatrix*Hill2Lvlh';
    end
end

CWMatrix(abs(CWMatrix)<1e-10) = 0;
end

