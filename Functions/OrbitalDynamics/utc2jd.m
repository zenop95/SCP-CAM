function jd = utc2jd(utc)
% utc2jd transforms a utc date into a julian date
%]
% INPUT: utc: utc date
%
% OUTPUT jd: [days] julian date

% Author: Zeno Pavanello 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
if isstring(utc) || isa(utc,"double")
    
elseif ischar(utc)
    utc(utc < '0' | utc > '9') = ' ';
    utc = sscanf(utc,'%f');
else
    error('utc input must be char or string')
end
y  = double(utc(3));
m  = double(utc(2)); 
d  = double(utc(1));
hh = double(utc(4));
mm = double(utc(5));
ss = double(utc(6));

hh = hh + mm/60 + ss/3600;
day = 1721013.5 + fix(275/9*m) + d + y*367 - ...
        fix(7/4*(y + fix((m+9)/12))); 
jd  = day + hh/24;
end