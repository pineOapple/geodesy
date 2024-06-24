function PLOT_POINT_SKYVIEW(azimuth, elevation)
    % Function to display a point on a Mollweide projection map
    % given the azimuth and elevation
    
    % Clear workspace and set up figure
    figure(2); clf(2);
    
    % Set up colormap
    mymap = colormap('jet'); 
    mymap(1,1:3) = 1;
    colormap(mymap);
    
    % Set up map axes
    axesm('MapProjection', 'mollweid', ...
          'MapLatLimit', [-90 90], ...
          'Gcolor', 'white', ...
          'GLineWidth', 2.0, ...
          'MLineLocation', [-135 -90 -45 0 45 90 135], ...
          'PLineLocation', 30);
    axis off;
    clim([0 500]);
    grid on;
    
    % Configure latitude and longitude labels
    plabel('LabelFormat', 'none');
    mlabel('on');
    parallel = 'equator';
    mlabel(parallel);
    mlabel('FontColor', 'white');
    mlabel('off');
    gridm('on');
    
    % Colorbar configuration
    c = colorbar('location', 'southoutside', ...
                 'box', 'on', ...
                 'color', [0 0 0]);
    c.Label.String = 'T_{sky} (K)';
    c.Limits = [0 500];
    c.Ticks = 0:50:500;
    c.FontSize = 12;
    
    % Add text annotations for time labels
    textm(-5, -135, '3 h', 'color', 'white', 'fontsize', 25);
    textm(-5, -90, '6 h', 'color', 'white', 'fontsize', 25);
    textm(-5, -45, '9 h', 'color', 'white', 'fontsize', 25);
    textm(-5, 0, '12 h', 'color', 'white', 'fontsize', 25);
    textm(-5, 45, '15 h', 'color', 'white', 'fontsize', 25);
    textm(-5, 90, '18 h', 'color', 'white', 'fontsize', 25);
    textm(-5, 135, '21 h', 'color', 'white', 'fontsize', 25);
    
    % Plot the point using the provided azimuth and elevation
    plotm(elevation, azimuth, 'r.', 'MarkerSize', 20);
    
    hold off;
end
