function [eps0, deleps, delpsi] = GDS_NUTATION(t)
% GDS_NUTATION Computes the mean obliquity of the ecliptic, nutation in longitude and obliquity
%
% DESCRIPTION:
%   This function calculates the mean obliquity of the ecliptic (epsilon_0),
%   the nutation in longitude (delta_psi), and the nutation in obliquity 
%   (delta_epsilon) for a given time since epoch J2000.0 in Julian centuries.
%
% USAGE:
%   [eps0, deleps, delpsi] = GDS_NUTATION(t)
%
% INPUT:
%   t - Time since epoch J2000.0 (in Julian centuries of 36525 days)
%
% OUTPUT:
%   eps0   - Mean obliquity of the ecliptic (in decimal degrees)
%   deleps - Nutation in obliquity (in decimal degrees)
%   delpsi - Nutation in longitude (in decimal degrees)
%
% NOTE:
%   The accuracy of the formulas is approximately one arcsecond.
%
% EXAMPLE:
%   t = 0.2; % Time since epoch J2000.0 in Julian centuries
%   [eps0, deleps, delpsi] = GDS_NUTATION(t);
%
% COPYRIGHT:
%   (c) Sneeuw/Zebhauser, IAPG, TU Munich
%
% AUTHOR:
%   Sneeuw/Zebhauser, IAPG, TU Munich
%
% DATE:
%   2024-06-20
%
%----------------------------------------------------------------------------

% Zeitdifferenz von J2000
d = t * 36525;

% Hilfsfunktionen
f1 = (125 - 0.05295 * d) / 180 * pi;
f2 = (200.9 + 1.97129 * d) / 180 * pi;

eps0 = (84381.448 - 46.815 * t) / 3600;
delpsi = (-0.0048 * sin(f1) - 0.0004 * sin(f2));
deleps = (0.0026 * cos(f1) + 0.0002 * cos(f2));

end
