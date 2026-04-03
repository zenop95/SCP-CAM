function mjd = utc2mjd(utc)
% utc2jd transforms a utc date into a modified julian date
%]
% INPUT: utc: utc date
%
% OUTPUT mjd: [days] modified julian date

% Author: Zeno Pavanello 2022
% E-mail: zpav176@aucklanduni.ac.nz
%-------------------------------------------------------------------------
mjd = (utc2jd(utc) - 2451545.5);
end