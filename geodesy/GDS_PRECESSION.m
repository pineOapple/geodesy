function [za, thetaa, zetaa] = GDS_PRECESSION(t)
% PREZWINK Computes the three precession angles z_A, theta_A, and zeta_A
%
% DESCRIPTION:
%   This function calculates the three precession angles (z_A, theta_A, and
%   zeta_A) for a given time since epoch J2000.0 in Julian centuries.
%
% USAGE:
%   [za, thetaa, zetaa] = prezwink(t)
%
% INPUT:
%   t - Time since epoch J2000.0 (in Julian centuries of 36525 days)
%
% OUTPUT:
%   za     - Precession angle z_A (in decimal degrees)
%   thetaa - Precession angle theta_A (in decimal degrees)
%   zetaa  - Precession angle zeta_A (in decimal degrees)
%
% NOTE:
%   The accuracy of the formulas is approximately one arcsecond.
%
% EXAMPLE:
%   t = 0.2; % Time since epoch J2000.0 in Julian centuries
%   [za, thetaa, zetaa] = prezwink(t);
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

zetaa  = (2306.2181 * t + 0.30188 * t.^2) / 3600;
za     = (2306.2181 * t + 1.09468 * t.^2) / 3600;
thetaa = (2004.3109 * t - 0.42665 * t.^2) / 3600;

end
