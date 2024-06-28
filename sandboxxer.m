clear
clc
close all

% Initialize the simulation with Neumayer Station's coordinates
phi = -70.6734; % latitude in degrees
lambda = -8.2741; % longitude in degrees
phi = 45; % latitude in degrees
lambda = 90; % longitude in degrees
phi = 0; % latitude in degrees
lambda = 89; % longitude in degrees
yr = 2024;
m = 3;
d = 21;
ut1 = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INERTIAL
figure;
set(gcf, 'Position',  [1720, -7, 1720, 1440]);
tiledlayout(2,3);
nexttile;
hold on;
view(30,30);
axis equal;

% [1] Make a sphere, plot it
[u, v, w] = ellipsoid(0,0,0,1,1,1,20);
surf(u, v, w, 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% [2] Create the GC and plot them
theta = linspace(-pi,pi,100);
x = cos(theta); y = sin(theta); z = zeros(size(theta));
gc_xy_i = [x' y' z']';
gc_xz_i = [x' z' y']';
gc_yz_i = [z' y' x']';
plot3(gc_xy_i(1,:), gc_xy_i(2,:), gc_xy_i(3,:), 'k', 'LineWidth', 1);
plot3(gc_xz_i(1,:), gc_xz_i(2,:), gc_xz_i(3,:), 'k', 'LineWidth', 1);
plot3(gc_yz_i(1,:), gc_yz_i(2,:), gc_yz_i(3,:), 'k', 'LineWidth', 1);

% [3] Create the GWC and plot it
theta = theta/2;
x = cos(theta); y = zeros(size(theta)); z = sin(theta);
gw_i = [x' y' z']';
plot3(gw_i(1,:), gw_i(2,:), gw_i(3,:), 'r', 'LineWidth', 2);

% [3] Create the Sunposition and plot it
[t, ~, ~] = GDS_JULIANC(yr, m, d, ut1);
[aapp,dapp,~] = GDS_SOLARPOS(t);
[x y z] = sph2cart(deg2rad(aapp),deg2rad(dapp),1);
pos = POINT_TO_VECTOR(x,y,z);
HOLD_PLOT_VEC3(pos, 'Sonnenposition: '+string(ut1)+'h');

% [5] Create the Ecliptic and plot it
gc_ec_i = rotx(23.5)*gc_xy_i; % Positive, cause we rotate the object!
plot3(gc_ec_i(1,:), gc_ec_i(2,:), gc_ec_i(3,:), 'k', 'LineWidth', 2);

% [6] Create the Ecliptic and plot it
[R, G] = GDS_INT_TO_LCL(lambda,phi,yr,m,d,0);
gw_c = G*gw_i;
h = plot3(gw_c(1,:), gw_c(2,:), gw_c(3,:), 'b', 'LineWidth', 2);

% [4] Create inertial transformation matrix and rotate GWC
for j = 1:1:1
    [R, G] = GDS_INT_TO_LCL(lambda,phi,yr,m,d,j);
    gw_c = G*gw_i;
    set(h, 'XData', gw_c(1, :), 'YData', gw_c(2, :), 'ZData', gw_c(3, :));
    pause(0.05);
end

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SKYPLOT
nexttile;
GDS_SKYPLOT(0, 0);
hold on; % Ensure all plots are on the same figure

colors = [
    "#f21821",  % Red
    "#f8631f",
    "#fa931a",
    "#ffc309",
    "#fff600",
    "#cdde25",
    "#8bc83b",
    "#04b99e",
    "#01aef3",
    "#5954a8",
    "#8f59a7",  % Purple
    "#bf168d"   % Magenta
];

% Array to store plot handles for legend
legend_handles = [];
legend_entries = {};

% Cell array of month names
month_names = {'January', 'February', 'March', 'April', 'May', 'June', ...
               'July', 'August', 'September', 'October', 'November', 'December'};

% [3] Create the Sunposition and plot it

for month = 1:1:6

    ut1 = 0;
    [R, G] = GDS_INT_TO_LCL(lambda, phi, yr, month, d, ut1);
    [t, ~, ~] = GDS_JULIANC(yr, month, d, ut1);
    [aapp, dapp, ~] = GDS_SOLARPOS(t);
    [x, y, z] = sph2cart(deg2rad(aapp), deg2rad(dapp), 1);
    vec = R * [x y z]';
    [az, el, ~] = cart2sph(vec(1), vec(2), vec(3));
    sunpos = [rad2deg(az) rad2deg(el)];
    
    for ut1 = 0:0.1:24
        [R, G] = GDS_INT_TO_LCL(lambda, phi, yr, month, d, ut1);
        [t, ~, ~] = GDS_JULIANC(yr, month, d, ut1);
        [aapp, dapp, ~] = GDS_SOLARPOS(t);
        [x, y, z] = sph2cart(deg2rad(aapp), deg2rad(dapp), 1);
        vec = R * [x y z]';
        [az, el, ~] = cart2sph(vec(1), vec(2), vec(3));
        sunpos(end+1, :) = [rad2deg(az) rad2deg(el)];
    end
    
    sunpos = sunpos(sunpos(:, 2) >= 0, :); % Filter out negative elevations
    
    % Plot and store handle
    hsky = GDS_SKYPLOT(sunpos(:, 1), sunpos(:, 2), '-', 3, colors(month, :), '');
    legend_handles(end+1) = hsky; % Collect handles
    legend_entries{end+1} = month_names{month}; % Collect legend entries
end

% Update legend at the end
legend(legend_handles, legend_entries);

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SKYPLOT
nexttile;
GDS_SKYPLOT(0, 0);
hold on; % Ensure all plots are on the same figure

colors = [
    "#f21821",  % Red
    "#f8631f",
    "#fa931a",
    "#ffc309",
    "#fff600",
    "#cdde25",
    "#8bc83b",
    "#04b99e",
    "#01aef3",
    "#5954a8",
    "#8f59a7",  % Purple
    "#bf168d"   % Magenta
];

% Array to store plot handles for legend
legend_handles = [];
legend_entries = {};

% Cell array of month names
month_names = {'January', 'February', 'March', 'April', 'May', 'June', ...
               'July', 'August', 'September', 'October', 'November', 'December'};

% [3] Create the Sunposition and plot it

for month = 7:1:12

    ut1 = 0;
    [R, G] = GDS_INT_TO_LCL(lambda, phi, yr, month, d, ut1);
    [t, ~, ~] = GDS_JULIANC(yr, month, d, ut1);
    [aapp, dapp, ~] = GDS_SOLARPOS(t);
    [x, y, z] = sph2cart(deg2rad(aapp), deg2rad(dapp), 1);
    vec = R * [x y z]';
    [az, el, ~] = cart2sph(vec(1), vec(2), vec(3));
    sunpos = [rad2deg(az) rad2deg(el)];
    
    for ut1 = 0:0.1:24
        [R, G] = GDS_INT_TO_LCL(lambda, phi, yr, month, d, ut1);
        [t, ~, ~] = GDS_JULIANC(yr, month, d, ut1);
        [aapp, dapp, ~] = GDS_SOLARPOS(t);
        [x, y, z] = sph2cart(deg2rad(aapp), deg2rad(dapp), 1);
        vec = R * [x y z]';
        [az, el, ~] = cart2sph(vec(1), vec(2), vec(3));
        sunpos(end+1, :) = [rad2deg(az) rad2deg(el)];
    end
    
    sunpos = sunpos(sunpos(:, 2) >= 0, :); % Filter out negative elevations
    
    % Plot and store handle
    hsky = GDS_SKYPLOT(sunpos(:, 1), sunpos(:, 2), '-', 3, colors(month, :), '');
    legend_handles(end+1) = hsky; % Collect handles
    legend_entries{end+1} = month_names{month}; % Collect legend entries
end

% Update legend at the end
legend(legend_handles, legend_entries);

hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL

nexttile;
hold on;
view(30,30);
axis equal;
% [1] Make a sphere, plot it
[u, v, w] = ellipsoid(0,0,0,1,1,1,20);
surf(u, v, w, 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% [2] Create the GC and plot them
theta = linspace(-pi,pi,100);
x = cos(theta); y = sin(theta); z = zeros(size(theta));
gc_xy_i = [x' y' z']';
gc_xz_i = [x' z' y']';
gc_yz_i = [z' y' x']';
plot3(gc_xy_i(1,:), gc_xy_i(2,:), gc_xy_i(3,:), 'k', 'LineWidth', 1);
plot3(gc_xz_i(1,:), gc_xz_i(2,:), gc_xz_i(3,:), 'k', 'LineWidth', 1);
plot3(gc_yz_i(1,:), gc_yz_i(2,:), gc_yz_i(3,:), 'k', 'LineWidth', 1);

% [3] Create the GWC and plot it
theta = theta/2;
x = cos(theta); y = zeros(size(theta)); z = sin(theta);
gw_i = [x' y' z']';
plot3(gw_i(1,:), gw_i(2,:), gw_i(3,:), 'r', 'LineWidth', 2);

% [3] Create the Spring Point and plot it
[t, ~, ~] = GDS_JULIANC(yr, m, d, ut1);
[aapp,dapp,~] = GDS_SOLARPOS(t);
[x, y, z] = sph2cart(deg2rad(aapp),deg2rad(dapp),1);
pos = POINT_TO_VECTOR(x,y,z);
HOLD_PLOT_VEC3(pos, 'Spring Point: '+string(ut1)+'h');

% [4] Create inertial transformation matrix and rotate GWC
[R, G] = GDS_INT_TO_LCL(lambda,phi,yr,m,d,ut1);
gw_i = G*gw_i;
plot3(gw_i(1,:), gw_i(2,:), gw_i(3,:), 'b', 'LineWidth', 2);

% [5] Create the Ecliptic and plot it
gc_ec_i = rotx(23.5)*gc_xy_i;
plot3(gc_ec_i(1,:), gc_ec_i(2,:), gc_ec_i(3,:), 'k', 'LineWidth', 2);

hold off;