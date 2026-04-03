function [x_ecef,dcm_EI] = eci2ecefSimp(mjd,x_eci)
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
ERA = 2*pi*jdFraction;
% Transformation matrix for earth rotation
R = eulang2rmat([ERA 0 0]', 'zyx')';

% Composition of CIP rotation and Earth rotation
dcm_EI = R;
% Position and velocity in ECEF coordinates
try
    x_ecef = [dcm_EI*x_eci(1:3,:); dcm_EI*x_eci(4:6,:)];
catch
    x_ecef = dcm_EI*x_eci(1:3,:);
end
end