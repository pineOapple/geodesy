clear
clc
close all

% Initialize the simulation with Neumayer Station's coordinates
phi = -70.6734; % latitude in degrees
lambda = -8.2741; % longitude in degrees
% phi = 0; % latitude in degrees
% lambda = 0; % longitude in degrees
season = 4;
alpha = [0 90 180 270];
delta = [0 23.5 0 -23.5];
cmonth = [3 6 9 12];


alpha = alpha(season);
delta = delta(season);
cmonth = cmonth(season);
cyear = 4224;
cmonth = 6;
cday = 21;
ut1 = 12;

figure;
set(gcf, 'Position',  [1720, -7, 1720, 1440]);
hold on;
view(30,30);
[u, v, w] = ellipsoid(0,0,0,1,1,1,20);
surf(u, v, w, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
axis equal;

% Define the GCs amd specific coordinates
theta = linspace(-pi, pi, 100);
x = cos(theta); y = sin(theta); z = zeros(size(theta));

gc_xy_i = [x' y' z']';
gc_xz_i = [x' z' y']';
gc_yz_i = [z' y' x']';

theta = theta/2;
x = cos(theta); y = zeros(size(theta)); z = sin(theta);
gw_i = [x' y' z']';

plot3(gc_xy_i(1,:), gc_xy_i(2,:), gc_xy_i(3,:), 'k', 'LineWidth', 1);
plot3(gc_xz_i(1,:), gc_xz_i(2,:), gc_xz_i(3,:), 'k', 'LineWidth', 1);
plot3(gc_yz_i(1,:), gc_yz_i(2,:), gc_yz_i(3,:), 'k', 'LineWidth', 1);
% for i=1:1:24
[t, ~, ~] = GDS_JULIANC(cyear, cmonth, cday, ut1);
[gast, ~] = GDS_JUL_TO_GAST(ut1, t);
[x,y,z] = sph2cart(deg2rad(gast*15), theta/2, 1);

plot3(gw_i(1,:), gw_i(2,:), gw_i(3,:), 'r', 'LineWidth', 2);
% end

% Plot stuff
pos = POINT_TO_VECTOR(1,0,0);
HOLD_PLOT_VEC3(pos, 'Spring Point: '+string(ut1)+'h');

% Simulate sun position over a day
sun_pos = zeros(24,4);
for i = 1:24
    sun_pos(i, 1:3) = INERTIAL_TO_LOCAL(lambda, phi, cyear, cmonth, cday, i, alpha, delta, 1);
    sun_pos(i, 4) = i;
end

sun_pos(i+1, :) = sun_pos(1, :);
sun_pos(i+1, 4) = 0;

plot3(sun_pos(:,1), sun_pos(:,2), sun_pos(:,3), 'b-', LineWidth=2);
text(sun_pos(:,1), sun_pos(:,2), sun_pos(:,3), string(sun_pos(:,4)));

hold off;
