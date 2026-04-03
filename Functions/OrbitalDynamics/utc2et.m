function et = utc2et(utc)
% utc2et transforms a utc date into a ephemeris time in seconds
%]
% INPUT: utc: utc date
%
% OUTPUT et: [s] ephemeris time

% Author: Zeno Pavanello 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
et = (utc2jd(utc) - 2451545)*86400;
end