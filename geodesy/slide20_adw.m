% MATLAB Script to visualize the position, orbit of Jupiter, and the ecliptic plane in the Solar System

% Constants for Jupiter's orbit (epoch J2000.0)
a = 5.2044;         % Semi-major axis in AU
e = 0.0489;         % Eccentricity
i = 1.303;          % Inclination in degrees
Omega = 100.464;    % Longitude of the ascending node in degrees
omega = 273.867;    % Argument of perihelion in degrees
M0 = 19.65053;      % Mean anomaly at epoch J2000.0 in degrees
n = 0.083056;       % Mean daily motion in degrees per day

% Define the time for which you want the ephemeris data
% This is the number of days since the standard epoch J2000.0 (2000-01-01 12:00:00)
date_vec = [2024, 07, 15, 0, 0, 0];
date_num = datenum(date_vec);
epoch_num = datenum([2000, 01, 01, 12, 0, 0]);
d = date_num - epoch_num;

% Calculate the current mean anomaly
M = M0 + n * d;

% Solve Kepler's equation to get the eccentric anomaly
E = M;  % Initial guess for E
for j = 1:1000
    E = M + e * sind(E);  % Update E using the iterative method
end

% Calculate the true anomaly
v = 2 * atan2d(sqrt(1 + e) * sind(E / 2), sqrt(1 - e) * cosd(E / 2));

% Calculate the heliocentric distance
r = a * (1 - e * cosd(E));

% Calculate the position in the orbital plane
x_orbit = r * cosd(v);
y_orbit = r * sind(v);

% Rotate to the ecliptic plane
x_ecl = x_orbit * (cosd(omega) * cosd(Omega) - sind(omega) * cosd(i) * sind(Omega)) - y_orbit * (sind(omega) * cosd(Omega) + cosd(omega) * cosd(i) * sind(Omega));
y_ecl = x_orbit * (cosd(omega) * sind(Omega) + sind(omega) * cosd(i) * cosd(Omega)) + y_orbit * (cosd(omega) * cosd(i) * cosd(Omega) - sind(omega) * sind(Omega));
z_ecl = x_orbit * sind(omega) * sind(i) + y_orbit * cosd(omega) * sind(i);

% Calculate the full orbit path
theta = linspace(0, 360, 1000);
r_orbit = a * (1 - e^2) ./ (1 + e * cosd(theta));
x_orbit_path = r_orbit .* cosd(theta);
y_orbit_path = r_orbit .* sind(theta);

% Rotate the orbit path to the ecliptic plane
x_ecl_path = x_orbit_path * (cosd(omega) * cosd(Omega) - sind(omega) * cosd(i) * sind(Omega)) - y_orbit_path * (sind(omega) * cosd(Omega) + cosd(omega) * cosd(i) * sind(Omega));
y_ecl_path = x_orbit_path * (cosd(omega) * sind(Omega) + sind(omega) * cosd(i) * cosd(Omega)) + y_orbit_path * (cosd(omega) * cosd(i) * cosd(Omega) - sind(omega) * sind(Omega));
z_ecl_path = x_orbit_path * sind(omega) * sind(i) + y_orbit_path * cosd(omega) * sind(i);

% Visualize the Solar System
figure;
hold on;
axis equal;
grid on;

% Plot the Sun
plot3(0, 0, 0, 'yo', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
text(1, 0, 0, '  Sun', 'FontSize', 12);

% Plot Jupiter's orbit
plot3(x_ecl_path, y_ecl_path, z_ecl_path, 'b-', 'LineWidth', 1.5);

% Plot Jupiter's current position
plot3(x_ecl, y_ecl, z_ecl, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(x_ecl, y_ecl, z_ecl, '  Jupiter', 'FontSize', 12);

% Plot the ecliptic plane
[X, Y] = meshgrid(linspace(-6, 6, 2), linspace(-6, 6, 2));
Z = zeros(size(X));
surf(X, Y, Z, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0, 1, 0]);
text(5, 0, 0, '  Ecliptic Plane', 'FontSize', 12, 'Color', 'g');

% Label the axes
xlabel('X (AU)');
ylabel('Y (AU)');
zlabel('Z (AU)');
title('Position and Orbit of Jupiter in the Solar System');

% Set the view
view(3);

hold off;
