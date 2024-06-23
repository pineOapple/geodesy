function [gast, gmst] = GDS_JUL_TO_GAST(ut1, t)
% JUL2GAST Computes the Greenwich Actual (true) Siderial Time (GAST)
%
% DESCRIPTION:
%   This function returns the Greenwich Actual Siderial Time (GAST) for a given
%   Universal Time (UT1) and time since epoch J2000.0 in Julian centuries.
%
% USAGE:
%   [gast, gmst] = jul2gast(ut1, t)
%
% INPUT:
%   ut1 - Universal Time in Greenwich (decimal hours)
%   t   - Time since epoch J2000.0 (in Julian centuries of 36525 days)
%
% OUTPUT:
%   gast - Greenwich Actual Siderial Time (decimal hours)
%   gmst - Greenwich Mean Siderial Time (decimal hours)
%
% NOTE:
%   Due to the accuracy level of nutation, the result will be accurate to
%   about 0.1 seconds (of time).
%
% EXAMPLE:
%   ut1 = 12.0; % Universal Time in decimal hours
%   t = 0.2;    % Time since epoch J2000.0 in Julian centuries
%   [gast, gmst] = jul2gast(ut1, t);
%
% COPYRIGHT:
%   (c) Sneeuw/Zebhauser, IAPG, TU Munich
%
% AUTHOR:
%   Sneeuw/Zebhauser, IAPG, TU Munich
%
% DATE:
%   1996-01-04
%
% REVISION HISTORY:
%   .Oct.2001, NS: Translation of original PREZWINK.M
%
% REMARKS:
%   This function uses the NUTATION subfunction to calculate nutation components.
%
%----------------------------------------------------------------------------

% Here we go
gmst = (ut1 * 3600 + 24110.54841 + 8640184.812866 * t + 0.093104 * t.^2 ...
        - 6.2e-6 * t.^3) / 3600;
[eps0, deleps, delpsi] = geodesy.GDS_NUTATION(t);
epsi = (eps0 + deleps) / 180 * pi;
eq = delpsi .* cos(epsi);
gast = gmst + eq / 15;
gast = rem(rem(gast, 24) + 24, 24);
gmst = rem(rem(gmst, 24) + 24, 24);
end
