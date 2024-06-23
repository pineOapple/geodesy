%% PREAMBEL
clear;
clc;
close all;

%% INIT
a = 1* 6378137; 
b = 1*6356752.3;
r = (a + b) / 2;
e = ELLIPSE_EXCENTRICITY(a, b);
E = ELLIPSE_EXCENTRICITY_LIN(a, e);
f = ELLIPSE_FLATTENING_1(a, b);
F = 1 / f;
phi = COORDINATE_TO_ANGLE(48, 46, 39) * 2 * pi / 360;
% phi = linspace(-pi / 2, pi / 2, 10);
lambda = phi;

%% PROCESS
[x_1, y_1] = ELLIPSE_X_Y(a, b, 100);
[x_2, y_2] = ELLIPSE_X_Y(a, a, 100);
[x_3, y_3] = ELLIPSE_X_Y(b, b, 100);
[x_m, y_m] = ELLIPSE_X_Y(E, E, 100);

%% POSTPROCESS
figure;
hold on;
grid on;
axis equal;
view(30, 30);

% Plot ellipses
plot3(zeros(size(x_1)), x_1, y_1, 'DisplayName', 'Ellipse a, b');
plot3(zeros(size(x_2)), x_2, y_2, 'DisplayName', 'Circle radius a');
plot3(zeros(size(x_3)), x_3, y_3, 'DisplayName', 'Circle radius b');
plot3(zeros(size(x_m)), x_m, y_m, 'DisplayName', 'Circle radius E');

% Labeled focal points
plot3(0, E, 0, 'ro', 'DisplayName', 'Focal Point 1');
plot3(0, -E, 0, 'ro', 'DisplayName', 'Focal Point 2');

% Center point
plot3(0, 0, 0, 'ro', 'DisplayName', 'Center');

% Calculate and plot radial lines
N = ELLIPSE_N(a, e, phi);
[x, y, z] = ELLIPSE_COORD(a, e, phi, 0);
x_end = N .* cos(phi);
y_end = N .* (1 - e^2) .* sin(phi);
x_start = zeros(size(phi));
y_start = -(N * (e^2)) .* sin(phi);

% GEOCETRIC 
psi = acos(x_end./sqrt(y_end.^2+x_end.^2));


for i = 1:length(phi)
    line([0 0], [x_start(i) x_end(i)], [y_start(i) y_end(i)], 'Color', 'black');
    line([0 0], [0 x_end(i)], [0 y_end(i)], 'Color', 'red');
    text(0, x_end(i), y_end(i), string(i)+ ' : '+ string(phi(i)*360/(2*pi)))
end

[u v w] = ellipsoid(0,0,0,a,a,b,50);
surf(u,v,w, 'EdgeColor', 'none', 'FaceAlpha', 0.1)
hold off;


'Exzentrizität: ' + string(e)
'Inverse Abplattung: ' + string(F)
'Radius Fokalkreis: ' + string(E/1000)
'Normalkrümmungsradius: ' + string(N/1000)
'z-Achseabschnitt (N): ' + string(min(y_start)/1000)
'Abstand Äquatorebene: ' + string(min(y_end)/1000)
'Abstand Rotationsachse: ' + string(min(x_end)/1000)
'Geodätische Breite: ' + string(phi.*360/(2*pi))
'Geozentrische Breite: ' + string(psi.*360/(2*pi))
[d1,m1,s1] = ANGLE_TO_COORD(rad2deg(phi))
[d2,m2,s2] = ANGLE_TO_COORD(rad2deg(psi))
'Geozentrische Abweichung: ' + sqrt(y_end.^2+x_end.^2)


grid on
%% Functions http://www.in-dubio-pro-geo.de/index.php?file=ellip/latit1


function e = ELLIPSE_EXCENTRICITY(a, b)
    e = sqrt((a^2 - b^2)) / a;
end

function E = ELLIPSE_EXCENTRICITY_LIN(a, e)
    E = a * e;
end

function [x, y] = ELLIPSE_X_Y(a, b, n)
    phi = linspace(0, 2 * pi, n);
    x = a * cos(phi);
    y = b * sin(phi);
end

function f = ELLIPSE_FLATTENING_1(a, b)
    f = (a - b) / b;
end

function phi = COORDINATE_TO_ANGLE(deg, min, sec)
    phi = ((sec / 60) + min) / 60 + deg;
end




