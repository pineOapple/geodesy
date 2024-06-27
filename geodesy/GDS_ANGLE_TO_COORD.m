function [deg, min, sec] = GDS_ANGLE_TO_COORD(phi)
% GDS_ANGLE_TO_COORD Converts a decimal degree angle to degrees, minutes, and seconds
%
% DESCRIPTION:
%   This function converts an angle in decimal degrees to degrees, minutes, 
%   and seconds.
%
% USAGE:
%   [deg, min, sec] = GDS_ANGLE_TO_COORD(phi)
%
% INPUT:
%   phi - Angle in decimal degrees (scalar or array)
%
% OUTPUT:
%   deg - Degrees (same size as phi)
%   min - Minutes (same size as phi)
%   sec - Seconds (same size as phi)
%
% EXAMPLE:
%   phi = 123.456; % Angle in decimal degrees
%   [deg, min, sec] = GDS_ANGLE_TO_COORD(phi);
%   % deg = 123, min = 27, sec = 21
%
% COPYRIGHT:
%   (c) 2024 Noel Ernsting Luz
%
% AUTHOR:
%   Noel Ernsting Luz
%
% DATE:
%   2024-06-20

    deg = floor(phi);
    min = floor((phi - deg) * 60);
    sec = floor(((phi - deg) * 60 - min) * 60);
end
