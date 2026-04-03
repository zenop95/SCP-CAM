function et = jd2et(jd)
% utc2et transforms a Julian date into a ephemeris time in seconds
%]
% INPUT: jd: Julian date
%
% OUTPUT et: [s] ephemeris time

% Author: Zeno Pavanello 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
et = (jd - 2451545)*86400;
end