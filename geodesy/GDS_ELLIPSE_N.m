function N = GDS_ELLIPSE_N(a, e, phi)
% GDS_ELLIPSE_N Computes the radius of curvature in the prime vertical
%
% DESCRIPTION:
%   This function calculates the radius of curvature in the prime vertical
%   (N) for a given geodetic latitude (phi), semi-major axis (a), and
%   eccentricity (e) of an ellipsoid.
%
% USAGE:
%   N = GDS_ELLIPSE_N(a, e, phi)
%
% INPUT:
%   a   - Semi-major axis (scalar)
%   e   - Eccentricity (scalar)
%   phi - Geodetic latitude (scalar or array in radians)
%
% OUTPUT:
%   N   - Radius of curvature in the prime vertical (same size as phi)
%
% EXAMPLE:
%   a = 6378137; % WGS-84 semi-major axis in meters
%   e = 0.08181919; % WGS-84 eccentricity
%   phi = deg2rad(45); % Convert 45 degrees to radians
%   N = GDS_ELLIPSE_N(a, e, phi);
%
% COPYRIGHT:
%   (c) 2024 Noel Ernsting Luz
%
% AUTHOR:
%   Noel Ernsting Luz
%
% DATE:
%   2024-06-20

    N = a ./ sqrt(1 - (e^2) * (sin(phi).^2));
end
