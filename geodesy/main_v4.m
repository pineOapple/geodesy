% MATLAB Skript zur Visualisierung der 6 keplerschen Bahnelemente

% Clear workspace and command window
clear;
clc;

% Definition der keplerschen Bahnelemente (Beispielwerte)
a = 1;           % Große Halbachse (in AE)
e = 0.1;         % Exzentrizität
i = 30;          % Inklination (in Grad)
Omega = 45;      % Länge des aufsteigenden Knotens (in Grad)
omega = 60;      % Argument des Perizentrums (in Grad)
M = 0;           % Mittlere Anomalie (in Grad)

% Konvertierung der Winkel in Bogenmaß
i_rad = deg2rad(i);
Omega_rad = deg2rad(Omega);
omega_rad = deg2rad(omega);
M_rad = deg2rad(M);

% Berechnung der exzentrischen Anomalie E aus der mittleren Anomalie M
E = M_rad;
for iter = 1:100  % Iteration zur Lösung der Kepler-Gleichung
    E = M_rad + e*sin(E);
end

% Berechnung der wahren Anomalie v
v = 2*atan(sqrt((1+e)/(1-e)) * tan(E/2));

% Berechnung der Position in der Orbitalebene
r = a * (1 - e^2) / (1 + e*cos(v));  % Radius
x_orb = r * cos(v);
y_orb = r * sin(v);
z_orb = 0;

% Rotationsmatrix für die Inklination
R_i = [1, 0, 0;
       0, cos(i_rad), -sin(i_rad);
       0, sin(i_rad), cos(i_rad)];

% Rotationsmatrix für die Länge des aufsteigenden Knotens
R_Omega = [cos(Omega_rad), -sin(Omega_rad), 0;
           sin(Omega_rad), cos(Omega_rad), 0;
           0, 0, 1];

% Rotationsmatrix für das Argument des Perizentrums
R_omega = [cos(omega_rad), -sin(omega_rad), 0;
           sin(omega_rad), cos(omega_rad), 0;
           0, 0, 1];

% Gesamte Rotationsmatrix
R = R_Omega * R_i * R_omega;

% Anwendung der Rotationsmatrix auf die Position
pos = R * [x_orb; y_orb; z_orb];

% Orbitalpunkte berechnen
theta = linspace(0, 2*pi, 100);
r = a * (1 - e^2) ./ (1 + e * cos(theta));
x_orb = r .* cos(theta);
y_orb = r .* sin(theta);
z_orb = zeros(size(theta));

% Anwendung der Rotationsmatrix auf die Orbitalpunkte
orbital_points = R * [x_orb; y_orb; z_orb];

% 3D-Plot der Umlaufbahn
figure;
plot3(orbital_points(1, :), orbital_points(2, :), orbital_points(3, :), 'b', 'LineWidth', 1.5);
hold on;
plot3(0, 0, 0, 'ro', 'MarkerFaceColor', 'r'); % Zentrum (Sonne)
plot3(pos(1), pos(2), pos(3), 'ko', 'MarkerFaceColor', 'k'); % Position des Körpers

% Plot der Ekliptik-Ebene
theta_ekliptik = linspace(0, 2*pi, 100);
x_ekliptik = cos(theta_ekliptik);
y_ekliptik = sin(theta_ekliptik);
z_ekliptik = zeros(size(theta_ekliptik));
fill3(x_ekliptik, y_ekliptik, z_ekliptik, 'g', 'FaceAlpha', 0.1);

% Darstellung der Bahnelemente

% Große Halbachse
plot3([0, a], [0, 0], [0, 0], 'k--', 'LineWidth', 1);
text(a/2, 0, 0, 'a', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Exzentrizität
focus_x = -a * e;
plot3([focus_x, 0], [0, 0], [0, 0], 'r--', 'LineWidth', 1);
text(focus_x/2, 0, 0, 'e', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

% Inklination
plot3([0, 0], [0, 0], [0, sin(i_rad)], 'm--', 'LineWidth', 1);
text(0, 0, sin(i_rad)/2, 'i', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Länge des aufsteigenden Knotens
plot3([0, cos(Omega_rad)], [0, sin(Omega_rad)], [0, 0], 'c--', 'LineWidth', 1);
text(cos(Omega_rad)/2, sin(Omega_rad)/2, 0, '\Omega', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

% Argument des Perizentrums
plot3([0, cos(omega_rad)], [0, sin(omega_rad)], [0, 0], 'y--', 'LineWidth', 1);
text(cos(omega_rad)/2, sin(omega_rad)/2, 0, '\omega', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

% Mittlere Anomalie (als wahre Anomalie dargestellt)
plot3([0, cos(v)], [0, sin(v)], [0, 0], 'g--', 'LineWidth', 1);
text(cos(v)/2, sin(v)/2, 0, 'M', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

% Plot-Eigenschaften
xlabel('X [AE]');
ylabel('Y [AE]');
zlabel('Z [AE]');
title('Visualisierung der keplerschen Bahnelemente');
grid on;
axis equal;
rotate3d on;
legend({'Orbitalbahn', 'Zentrum (Sonne)', 'Körper', 'Ekliptik', 'Große Halbachse (a)', 'Exzentrizität (e)', 'Inklination (i)', 'Länge des aufsteigenden Knotens (\Omega)', 'Argument des Perizentrums (\omega)', 'Mittlere Anomalie (M)'}, 'Location', 'best');
hold off;
