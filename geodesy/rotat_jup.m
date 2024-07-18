% MATLAB Skript zur Visualisierung von alpha0, delta0 und W

% Clear workspace and command window
clear;
clc;

% Parameterdefinitionen (Beispielwerte für einen Himmelskörper)
alpha0 = 257.31;  % Rektaszension des Nordpols in Grad
delta0 = -15.18;  % Deklination des Nordpols in Grad
W = 203.81;       % Längengrad des Primärmeridians in Grad

% Konvertierung der Winkel in Bogenmaß für Plotfunktionen
alpha0_rad = deg2rad(alpha0);
delta0_rad = deg2rad(delta0);
W_rad = deg2rad(W);

% Erstellen einer Kugel zur Darstellung des Himmelskörpers
[x, y, z] = sphere(50);

% Figur zur Darstellung des Himmelskörpers
figure;
surf(x, y, z, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
hold on;
axis equal;
colormap gray;
title('Visualisierung von alpha0, delta0 und W');

% Plot der Rotationsachse
quiver3(0, 0, 0, cos(delta0_rad) * cos(alpha0_rad), cos(delta0_rad) * sin(alpha0_rad), sin(delta0_rad), 1, 'r', 'LineWidth', 2);
text(cos(delta0_rad) * cos(alpha0_rad), cos(delta0_rad) * sin(alpha0_rad), sin(delta0_rad), '\alpha_0, \delta_0', 'FontSize', 12, 'Color', 'r');

% Plot des Primärmeridians
quiver3(0, 0, 0, cos(W_rad), sin(W_rad), 0, 1, 'b', 'LineWidth', 2);
text(cos(W_rad), sin(W_rad), 0, 'W', 'FontSize', 12, 'Color', 'b');

% Beschriftung der Achsen
xlabel('X');
ylabel('Y');
zlabel('Z');

% Großkreise für alpha (Längengrade)
alpha_values = deg2rad([0, 30, 60, 90, -30, -60, -90]);
theta = linspace(0, 2*pi, 100);
for alpha = alpha_values
    x_gc = cos(alpha) * cos(theta);
    y_gc = cos(alpha) * sin(theta);
    z_gc = sin(alpha) * ones(size(theta));
    plot3(x_gc, y_gc, z_gc, 'k--', 'LineWidth', 1);
end

% Großkreise für delta (Breitengrade)
delta_values = deg2rad([0, 30, 60, 90, -30, -60, -90]);
for delta = delta_values
    x_gc = cos(theta) .* cos(delta);
    y_gc = sin(theta) .* cos(delta);
    z_gc = sin(delta) * ones(size(theta));
    plot3(x_gc, y_gc, z_gc, 'k--', 'LineWidth', 1);
end

% Hinzufügen einer Legende
legend({'Himmelskörper', 'Rotationsachse', 'Primärmeridian', 'Längengrade', 'Breitengrade'}, 'Location', 'northeast');

% Achsenlimits festlegen
xlim([-1 1]);
ylim([-1 1]);
zlim([-1 1]);

% Grid einschalten
grid on;

% Ansicht einstellen
view(3);
rotate3d on;

hold off;
