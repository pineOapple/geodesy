function [x, y, z] = GDS_ELLIPSOID_COORD(a, e, phi, lambda)
% GDS_ELLIPSOID_COORD Converts geodetic coordinates to Cartesian coordinates
%
% DESCRIPTION:
%   This function converts geodetic coordinates (latitude, longitude) to 
%   Cartesian coordinates (x, y, z) for a given ellipsoid defined by the 
%   semi-major axis (a) and eccentricity (e).
%
% USAGE:
%   [x, y, z] = GDS_ELLIPSOID_COORD(a, e, phi, lambda)
%
% INPUT:
%   a      - Semi-major axis (scalar)
%   e      - Eccentricity (scalar)
%   phi    - Geodetic latitude (scalar or array in radians)
%   lambda - Geodetic longitude (scalar or array in radians)
%
% OUTPUT:
%   x - Cartesian coordinate X (same size as phi and lambda)
%   y - Cartesian coordinate Y (same size as phi and lambda)
%   z - Cartesian coordinate Z (same size as phi and lambda)
%
% EXAMPLE:
%   a = 6378137; % WGS-84 semi-major axis in meters
%   e = 0.08181919; % WGS-84 eccentricity
%   phi = deg2rad(45); % Convert 45 degrees latitude to radians
%   lambda = deg2rad(90); % Convert 90 degrees longitude to radians
%   [x, y, z] = GDS_ELLIPSOID_COORD(a, e, phi, lambda);
%
% COPYRIGHT:
%   (c) 2024 Noel Ernsting Luz
%
% AUTHOR:
%   Noel Ernsting Luz
%
% DATE:
%   2024-06-20

    N = GDS_ELLIPSE_N(a, e, phi); % Use the previously defined function
    x = N .* cos(phi) .* cos(lambda);
    y = N .* cos(phi) .* sin(lambda);
    z = N .* (1 - e^2) .* sin(phi);
end
