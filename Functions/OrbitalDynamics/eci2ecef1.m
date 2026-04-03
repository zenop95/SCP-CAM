function [x_ecef,dcm_EI] = eci2ecef(mjd,x_eci)
% Provides the Direction Cosine Matrix to rotate the ECI frame to the ECEF
% one. This function is based on the IAU-2000/2006 (IERS Technical Note 3).
% The computation neglects the Polar motion , the adjustment to the
% location of the Celestial Intermediate Pole, the difference between UTC
% and UT, assuming them as zero.
%
% INPUT:
%       mjd: the date/time in Modified Julian Date convention
% Parameters:
%       lunSolNut: numerical coefficients for lunar/solar nutation.
%       plaNut: numerical coefficients for planetary nutation.
%       polDev: numerical coefficients for polinomial evaluation in eq. 5.16.
%       sCIP: empirical coefficients used in the computation of S.
%       xCIP: empirical coefficients used in the computation of X.
%       yCIP: empirical coefficients used in the computation of Y.
% OUTPUT:
%       dcm: the direction cosine matrix
%
% References:
%   [1] Gerard Petit and Brian Luzum.,
%       "IERS Technical Note No. 36, IERS Conventions (2010).",
%       International Earth Rotation and Reference Systems Service (2010);
%       cap. 5, pages 43 - 73.
%   [2] FTP server for the numeric parameters:
%      	ftp://tai.bipm.org/iers/conv2010/chapter5/

% Authors: Zeno Pavanello, Davide Vertuani 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

%% CIO Based Transformation
% See IERS technical note 36 (table [5.2c])
% Elapsed Julian days since J2000
jdElapsed = mjd - 51544.5;

% Julian day fraction
jdFraction = mod(jdElapsed, 1);

%% Earth rotation
% Earth rotation angle (eq. [5.15])
ERA = mod(2*pi*(jdFraction + 0.7790572732640 + 0.00273781191135448*jdElapsed),2*pi);
% Transformation matrix for earth rotation
R = eulang2rmat([ERA 0 0]', 'zyx')';

%% Time and time powers
% Julian date for terrestrial time: 32.184 sec before jd
jdTerTime = mjd + 32.184/(24*3600);

% Number of Julian centuries since J2000 for terrestrial time.
tTime = (jdTerTime - 51544.5)/36525;

% Powers of the number of julian centuries
tTime2 = tTime^2;
tTime3 = tTime^3;
tTime4 = tTime^4;
tTime5 = tTime^5;

%% Celestial Motion of the CIP
load('eciEcef.mat','lunSolNut','plaNut','polDev','sCIP','xCIP','yCIP');
% Lunisolar nutation (eq. [5.43])
F       = zeros(14, 1);
F(1:5)  = lunSolNut(:, :)*[1 tTime tTime2 tTime3 tTime4]';
F       = deg2rad(F./3600)';

% Planetary nutation (eq. [5.44])
F(6:14) = plaNut(:,:)*[1 tTime]';

% In pa / F14 time is squared (eq. [5.44])
F(14) = F(14) * tTime;

% Vector arrangement
nutation = mod(F, 2*pi);
nutation = nutation';

% X-series (see eq. [5.16])
% tVec (tector with time powers) follows the separation in the provided
% coefficient table
tVecX              = zeros(length(xCIP),1);
tVecX(1:1306,:)    = ones(1306,1);
tVecX(1307:1559,:) = tTime  .* ones(253, 1);
tVecX(1560:1595,:) = tTime2 .* ones(36, 1);
tVecX(1596:1599,:) = tTime3 .* ones(4, 1);
tVecX(1600,:)      = tTime4;

argX = xCIP(:, 4:17) * nutation;
xSerie = sum((xCIP(:,2).*sin(argX(:,1)) + xCIP(:,3).*cos(argX(:,1))).*tVecX(:,1));

% Y-series (see eq. [5.16])
tVecY              = zeros(length(yCIP),1);
tVecY(1:962,:)     = ones(962,1);
tVecY(963:1239,:)  = tTime  .* ones(277, 1);
tVecY(1240:1269,:) = tTime2 .* ones(30, 1);
tVecY(1270:1274,:) = tTime3 .* ones(5, 1);
tVecY(1275,:)      = tTime4;

argY = yCIP(:,4:17) * nutation;
ySerie = sum((yCIP(:,2).*sin(argY(:,1)) + yCIP(:,3).*cos(argY(:,1))).*tVecY(:,1));

% S-series
tVecS          = zeros(length(sCIP),1);
tVecS(1:33,:)  = ones(33,1);
tVecS(34:36,:) = tTime  .* ones(3, 1);
tVecS(37:61,:) = tTime2 .* ones(25, 1);
tVecS(62:65,:) = tTime3 .* ones(4, 1);
tVecS(66,:)    = tTime4;

argS = sCIP(:,4:11) * [nutation(1:5,:); nutation(7:8,:); nutation(14,:)];
sSerie = sum((sCIP(:,2).*sin(argS(:,1)) + sCIP(:,3).*cos(argS(:,1))).*tVecS(:,1));

% Polynomial part of X and Y (eq. [5.16])
polinomial = polDev * [1 tTime tTime2 tTime3 tTime4 tTime5]';

% Adding the elements:
x = xSerie + polinomial(1);
y = ySerie + polinomial(2);
s = sSerie + polinomial(3);

% Microarcseconds to radians conversion
x = deg2rad(x*1e-6/3600);
y = deg2rad(y*1e-6/3600);
s = deg2rad(s*1e-6/3600) - x.*y/2;

% Coordinates of the CIP
E = atan2(y,x);
d = atan(sqrt((x.^2+y.^2)./(1-x.^2-y.^2)));
% Transformation matrix for celestial motion of the CIP
Q = eulang2rmat([E, d, -E-s]', 'zyz')';
% Composition of CIP rotation and Earth rotation
dcm_EI = R;%*Q;
% Position and velocity in ECEF coordinates
try
    x_ecef = [dcm_EI*x_eci(1:3,:); dcm_EI*x_eci(4:6,:)];
catch
    x_ecef = dcm_EI*x_eci(1:3,:);
end
end