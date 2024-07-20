clear
clc
close all

% Initialize the simulation with Neumayer Station's coordinates
phi     = -70.6734; % latitude in degrees
lambda  = -8.2741; % longitude in degrees
phi     = 0; % latitude in degrees
lambda  = 0; % longitude in degrees
season  = 4;
alpha   = [0 90 180 270];
delta   = [0 23.5 0 -23.5];
cmonth  = [3 6 9 12];
alpha   = alpha(season);
delta   = delta(season);
cmonth  = cmonth(season);
cyear   = 2024;
cmonth  = 6;
cday    = 21;
ut1     = 12;

% GAST, NUTATION, PRECESSION, JULICANCENTURIES
[t, ~, ~] = GDS_JULIANC(cyear, cmonth, cday, ut1);
[gast, ~] = GDS_JUL_TO_GAST(ut1, t);
[eps0, deleps, delpsi] = GDS_NUTATION(t);
[za, thetaa, zetaa] = GDS_PRECESSION(t);

% POSITION OF OBJECT
[r_i(1),r_i(2),r_i(3)] = sph2cart(deg2rad(alpha), deg2rad(delta), 1);

% TRANSFORMATION
N = rotx(eps0+deleps)*rotz(delpsi)*rotx(-eps0);
P = rotz(za)*roty(-thetaa)*rotz(zetaa);

r_1 = P*r_i';
r_2 = N*P*r_i';
r_3 = rotz(-(360*gast)/24)*N*P*r_i';
r_4 = rotz(-lambda)*rotz(-(360*gast)/24)*N*P*r_i';
r_5 = roty(-90+phi)*rotz(-lambda)*rotz(-(360*gast)/24)*N*P*r_i';

% Setup tiled layout for plotting
tiledlayout(2,3); % Define one row, two columns
set(gcf, 'Position',  [1720, -7, 1720, 1440])

% [ 1 1 ] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile;
hold on;

view(30,30);
[u, v, w] = ellipsoid(0,0,0,1,1,1,20);
surf(u, v, w, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
axis equal;

% Define the GCs
theta = linspace(0, 2*pi, 100);
x = cos(theta); y = sin(theta); z = zeros(size(theta));
plot3(x, y, z, 'k', 'LineWidth', 1);
plot3(x, z, y, 'k', 'LineWidth', 1);
plot3(z, y, x, 'k', 'LineWidth', 1);

theta = linspace(-pi/2, pi/2, 100);
% for i=1:1:24
[t, ~, ~] = GDS_JULIANC(cyear, cmonth, cday, ut1);
[gast, ~] = GDS_JUL_TO_GAST(ut1, t);
[x,y,z] = sph2cart(deg2rad(gast*15), theta, 1);
plot3(x, y, z, 'r', 'LineWidth', 1);
% end

% Plot stuff
pos = POINT_TO_VECTOR(1.4,0,0);
HOLD_PLOT_VEC3(pos, 'Spring Point: '+string(ut1)+'h');
HOLD_PLOT_VEC3(r_1,'Transformation: P');
hold off;

% [ 1 2 ] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile;
hold on;

view(30,30);
[u, v, w] = ellipsoid(0,0,0,1,1,1,20);
surf(u, v, w, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
axis equal;

% Define the GCs
theta = linspace(0, 2*pi, 100);
x = cos(theta); y = sin(theta); z = zeros(size(theta));
plot3(x, y, z, 'k', 'LineWidth', 1);
plot3(x, z, y, 'k', 'LineWidth', 1);
plot3(z, y, x, 'k', 'LineWidth', 1);

theta = linspace(-pi/2, pi/2, 100);
[t, ~, ~] = GDS_JULIANC(cyear, cmonth, cday, ut1);
[gast, ~] = GDS_JUL_TO_GAST(ut1, t);
[x,y,z] = sph2cart(deg2rad(gast*15), theta, 1);
plot3(x, y, z, 'r', 'LineWidth', 1);

% Plot stuff
pos = POINT_TO_VECTOR(1.4,0,0);
HOLD_PLOT_VEC3(pos, 'Spring Point: '+string(ut1)+'h');
HOLD_PLOT_VEC3(r_2,'Transformation: P+N');

hold off;

% [ 1 3 ] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile;
hold on;

view(30,30);
[u, v, w] = ellipsoid(0,0,0,1,1,1,20);
surf(u, v, w, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
axis equal;

% Define the GCs
theta = linspace(0, 2*pi, 100);
x = cos(theta); y = sin(theta); z = zeros(size(theta));
plot3(x, y, z, 'k', 'LineWidth', 1);
plot3(x, z, y, 'k', 'LineWidth', 1);
plot3(z, y, x, 'k', 'LineWidth', 1);
theta = linspace(-pi/2, pi/2, 100);
[t, ~, ~] = GDS_JULIANC(cyear, cmonth, cday, ut1);
[gast, ~] = GDS_JUL_TO_GAST(ut1, t);
[x,y,z] = sph2cart(deg2rad(gast*15), theta, 1);
plot3(x, y, z, 'r', 'LineWidth', 1);


% Plot stuff
pos = POINT_TO_VECTOR(1.4,0,0);
HOLD_PLOT_VEC3(pos, 'Spring Point: '+string(ut1)+'h');
HOLD_PLOT_VEC3(r_3,'Transformation: P+N+GAST');

hold off;

% [ 2 1 ] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile;
hold on;

view(30,30);
[u, v, w] = ellipsoid(0,0,0,1,1,1,20);
surf(u, v, w, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
axis equal;

% Define the GCs
theta = linspace(0, 2*pi, 100);
x = cos(theta); y = sin(theta); z = zeros(size(theta));
plot3(x, y, z, 'k', 'LineWidth', 1);
plot3(x, z, y, 'k', 'LineWidth', 1);
plot3(z, y, x, 'k', 'LineWidth', 1);
theta = linspace(-pi/2, pi/2, 100);
[t, ~, ~] = gds.GDS_JULIANC(cyear, cmonth, cday, ut1);
[gast, ~] = gds.GDS_JUL_TO_GAST(ut1, t);
[x,y,z] = sph2cart(deg2rad(gast*15), theta, 1);
plot3(x, y, z, 'r', 'LineWidth', 1);

% Plot stuff
pos = POINT_TO_VECTOR(1.4,0,0);
HOLD_PLOT_VEC3(pos, 'Spring Point: '+string(ut1)+'h');
HOLD_PLOT_VEC3(r_4,'Transformation: P+N+GAST+Lambda');

hold off;

% [ 2 2 ] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile;
hold on;

view(30,30);
[u, v, w] = ellipsoid(0,0,0,1,1,1,20);
surf(u, v, w, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
axis equal;

% Define the GCs
theta = linspace(0, 2*pi, 100);
x = cos(theta); y = sin(theta); z = zeros(size(theta));
plot3(x, y, z, 'k', 'LineWidth', 1);
plot3(x, z, y, 'k', 'LineWidth', 1);
plot3(z, y, x, 'k', 'LineWidth', 1);
theta = linspace(-pi/2, pi/2, 100);
[t, ~, ~] = gds.GDS_JULIANC(cyear, cmonth, cday, ut1);
[gast, ~] = gds.GDS_JUL_TO_GAST(ut1, t);
[x,y,z] = sph2cart(deg2rad(gast*15), theta, 1);
plot3(x, y, z, 'r', 'LineWidth', 1);

% Plot stuff
pos = POINT_TO_VECTOR(1.4,0,0);
HOLD_PLOT_VEC3(pos, 'Spring Point: '+string(ut1)+'h');
HOLD_PLOT_VEC3(r_5,'Transformation: P+N+GAST+Lambda+Phi');

hold off;

% [ 2 3 ] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Local results
nexttile;
hold on;
view(30,30);
axis equal;

[u, v, w] = ellipsoid(0,0,0,1,1,1,20);
surf(u, v, w, 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Define the GCs
theta = linspace(0, 2*pi, 100);
x = cos(theta); y = sin(theta); z = zeros(size(theta));
plot3(x, y, z, 'k', 'LineWidth', 1);
plot3(x, z, y, 'k', 'LineWidth', 1);
plot3(z, y, x, 'k', 'LineWidth', 1);
theta = linspace(-pi/2, pi/2, 100);
[t, ~, ~] = gds.GDS_JULIANC(cyear, cmonth, cday, ut1);
[gast, ~] = gds.GDS_JUL_TO_GAST(ut1, t);
[x,y,z] = sph2cart(deg2rad(gast*15), theta, 1);
plot3(x, y, z, 'r', 'LineWidth', 1);

% Simulate sun position over a day
sun_pos = zeros(24,4);
for k = 1:12
    for i = 1:24
        sun_pos(i, 1:3) = gds.INERTIAL_TO_LOCAL(lambda, phi, cyear, cmonth, cday, i, alpha, delta, 1);
        sun_pos(i, 4) = i;
    end
end
sun_pos(i+1, :) = sun_pos(1, :);
sun_pos(i+1, 4) = 0;

plot3(sun_pos(:,1), sun_pos(:,2), sun_pos(:,3), 'b-', LineWidth=2);
text(sun_pos(:,1), sun_pos(:,2), sun_pos(:,3), string(sun_pos(:,4)));

pos = gds.INERTIAL_TO_LOCAL(lambda, phi, cyear, cmonth, cday, ut1, 0, 0, 1.4);
HOLD_PLOT_VEC3(pos, 'Spring Point: '+string(ut1)+'h');

hold off;
