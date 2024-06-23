function hsky = GDS_SKYPLOT(azim, elev, line_style, line_width, line_color)
% GDS_SKYPLOT Polar coordinate plot using azimuth and elevation data
%
% DESCRIPTION:
%   This function creates a polar plot of the azimuth (AZIM) versus the elevation 
%   (ELEV). Negative elevation is allowed, and azimuth is counted clockwise from 
%   the North.
%
% USAGE:
%   hsky = GDS_SKYPLOT(azim, elev)
%   hsky = GDS_SKYPLOT(azim, elev, line_style, line_width, line_color)
%
% INPUT:
%   azim       - Azimuth angles in degrees (numeric array)
%   elev       - Elevation angles in degrees (numeric array)
%   line_style - Line style for the plot (optional, default: '-')
%   line_width - Line width for the plot (optional, default: 1.5)
%   line_color - Line color for the plot (optional, default: 'b')
%
% OUTPUT:
%   hsky - Handle to the plot object
%
% NOTE:
%   This function makes use of MATLAB's polar plotting capabilities. The default 
%   values for line style, line width, and line color are '-', 1.5, and 'b', 
%   respectively.
%
% EXAMPLE:
%   azim = 0:360; % Azimuth from 0 to 360 degrees
%   elev = 30 * ones(size(azim)); % Elevation at 30 degrees
%   hsky = GDS_SKYPLOT(azim, elev, '--', 2, 'r'); % Red dashed line
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

    % Default values
    if nargin < 3
        line_style = '-';
    end
    if nargin < 4
        line_width = 1.5;
    end
    if nargin < 5
        line_color = 'b';
    end

    % Checks and validations
    if ischar(azim) || ischar(elev)
        error('Input arguments must be numeric.');
    end
    if any(size(azim) ~= size(elev))
        error('AZIM and ELEV must be the same size.');
    end
    if any(elev > 90) || any(elev < -90)
        error('ELEV must be within [-90; 90] degrees.');
    end

    % Get hold state
    cax = newplot;
    next = lower(get(cax, 'NextPlot'));
    hold_state = ishold;

    % Get x-axis text color for grid lines
    tc = get(cax, 'xcolor');

    % Only draw grids if hold is off
    if ~hold_state

        % Make a radial grid
        hold on;
        % Check radial limits and ticks
        zenmax = max(90 - elev(:));
        zenmax = 15 * ceil(zenmax / 15);
        elmax  = 90;
        % Define a circle
        az = 0:pi/50:2*pi;
        xunit = sin(az);
        yunit = cos(az);

        % Draw solid circles each 30 degrees, with the horizon (elev=0) thicker
        for i = [30 60]
            plot(xunit * i, yunit * i, '-', 'color', tc, 'linewidth', 1);
        end
        i = 90; plot(xunit * i, yunit * i, '-', 'color', tc, 'linewidth', 2);
        for i = [15:30:75 105:15:zenmax]
            plot(xunit * i, yunit * i, ':', 'color', tc, 'linewidth', 1);
        end
        for i = 30:30:zenmax
            text(0, i, [' ' num2str(90 - i)], 'verticalalignment', 'bottom');
        end

        % Plot spokes
        az = (1:6) * 2 * pi / 12; % Define half circle
        caz = cos(az); 
        saz = sin(az);
        ca = [-caz; caz];
        sa = [-saz; saz];
        plot(elmax * ca, elmax * sa, '-', 'color', tc, 'linewidth', 1);
        if zenmax > elmax
            plot(zenmax * ca, zenmax * sa, ':', 'color', tc, 'linewidth', 1);
        end

        % Annotate spokes in degrees
        rt = 1.1 * elmax;
        for i = 1:length(az)
            loc1 = int2str(i * 30);
            if i == length(az)
                loc2 = int2str(0);
            else
                loc2 = int2str(180 + i * 30);
            end
            text(rt * saz(i), rt * caz(i), loc1, 'horizontalalignment', 'center');
            text(-rt * saz(i), -rt * caz(i), loc2, 'horizontalalignment', 'center');
        end

        % Brush up axes
        view(0, 90);
        axis(max(zenmax, elmax) * [-1 1 -1.1 1.1]);
        set(cax, 'position', [.05 .05 .9 .9])
    end

    % Transform data to Cartesian coordinates
    yy = (90 - elev) .* cos(azim / 180 * pi);
    xx = (90 - elev) .* sin(azim / 180 * pi);

    % Plot data on top of grid 
    q = plot(xx, yy, line_style, 'LineWidth', line_width, 'Color', line_color);
    if nargout > 0, hsky = q; end

    if ~hold_state
        axis('equal')
        axis('off')
    end

    % Reset hold state
    if ~hold_state, set(cax, 'NextPlot', next); end
end
