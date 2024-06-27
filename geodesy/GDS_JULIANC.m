function [t, mjd, jd] = GDS_JULIANC(y, m, d, ut1)
% JULIANJH Computes the time since epoch J2000 in Julian centuries
%
% DESCRIPTION:
%   This function calculates the time since the epoch J2000.0 in Julian 
%   centuries, the Modified Julian Date (MJD), and the Julian Date (JD).
%
% USAGE:
%   [t, mjd, jd] = julianjh(y, m, d, ut1)
%
% INPUT:
%   y   - Year (full year, e.g., 2024)
%   m   - Month (1-12)
%   d   - Day (1-31)
%   ut1 - Universal Time in Greenwich (decimal hours)
%
% OUTPUT:
%   t   - Time since epoch J2000.0 in Julian centuries (36525 days)
%   mjd - Modified Julian Date (days)
%   jd  - Julian Date (days)
%
% NOTE:
%   The conversion is exact.
%
% EXAMPLE:
%   y = 2024; % Year
%   m = 6;    % Month
%   d = 20;   % Day
%   ut1 = 12.0; % Universal Time in decimal hours
%   [t, mjd, jd] = julianjh(y, m, d, ut1);
%
% COPYRIGHT:
%   (c) Sneeuw/Zebhauser, IAPG, TU Munich
%
% AUTHOR:
%   Sneeuw/Zebhauser, IAPG, TU Munich
%
% DATE:
%   
%
% REVISION HISTORY:
%   - 980121NS: Adjusted y/m/d checks. Allowed vector input for y/m/d.
%   - 980122NS: Added MJD output.
%
%----------------------------------------------------------------------------

if any(m(:) > 12 | m(:) < 1)
    error('Monat zwischen 1 und 12')
end
if any(d(:) > 31 | d(:) < 1)
    error('Tag zwischen 1 und 31')
end
if any(rem(y(:), 1) ~= 0) || any(rem(m(:), 1) ~= 0) || any(rem(d(:), 1) ~= 0)
    error('Jahr, Monat, Tag muessen ganzzahlig eingegeben werden')
end

% Zwischenvariable JD ist das Julianische Datum
jd = 367 * y - floor(7 * (y + floor((m + 9) / 12)) / 4);
jd = jd + floor(275 * m / 9) + d + 1721014 + ut1 / 24 - 0.5;
t = (jd - 2451545) / 36525;
mjd = jd - 2400000.5; % modifiziertes jd

end
