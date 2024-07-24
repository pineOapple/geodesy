clear
clc
close all

% Einlesen files
filename = 'data/eop_iau2000_1846_now.txt';

% Laden der Daten
data = readtable(filename, 'Format','%f%f%f%f%f%f%f%f%f%f%f');
year = data{:, 1};
xP = data{:, 2};
xP_err = data{:, 3};
yP = data{:, 4};
yP_err = data{:, 5};

t = year;

% 1. Visualisierung der Zeitreihen
figure;
set(gcf, 'Position',  [500, 100, 2400, 600]);

colormap("hot")
subplot(2,1,1);
plot(t, xP);
title('Zeitreihe Polkoordinate x_P');
xlabel('Zeit (Jahre)');
ylabel('xP (arcsec)');

subplot(2,1,2);
plot(t, yP);
title('Zeitreihe Polkoordinate y_P');
xlabel('Zeit (Jahre)');
ylabel('yP (arcsec)');

saveas(gcf, 'plots/zeitreihen.png');

% 4. Spektralanalyse
Fs = 1 / mean(diff(t)); % Abtastfrequenz
xP_detrended = detrend(xP); % Trend entfernen
yP_detrended = detrend(yP); % Trend entfernen
n = length(t);
f = (0:n-1)*(Fs/n);
xP_fft = fft(xP_detrended);
yP_fft = fft(yP_detrended);
xP_power = abs(xP_fft).^2/n;
yP_power = abs(yP_fft).^2/n;
[~, xP_max_idx] = max(xP_power(2:end/2));
[~, yP_max_idx] = max(yP_power(2:end/2));
xP_main_period = 1/f(xP_max_idx + 1); % +1 wegen Index-Offset
yP_main_period = 1/f(yP_max_idx + 1);
period = 1 ./ f;

figure;
subplot(2,1,1);
plot(period(2:end/2), xP_power(2:end/2));
set(gca, 'XScale', 'log'); % Logarithmische Skala für bessere Darstellung
title('Leistungsspektrum Polkoordinate x_P');
xlabel('Periode (Jahre)');
ylabel('Leistung');

subplot(2,1,2);
plot(period(2:end/2), yP_power(2:end/2));
set(gca, 'XScale', 'log'); % Logarithmische Skala für bessere Darstellung
title('Leistungsspektrum Polkoordinate y_P');
xlabel('Periode (Jahre)');
ylabel('Leistung');
saveas(gcf, 'plots/leistungsspektrum.png');

% 5. Säkulare Polbewegung
p_xP = polyfit(t, xP, 1);
xP_trend = polyval(p_xP, t);
dxP_dt = p_xP(1); % Änderung pro Jahr

p_yP = polyfit(t, yP, 1);
yP_trend = polyval(p_yP, t);
dyP_dt = p_yP(1); % Änderung pro Jahr

figure;
set(gcf, 'Position',  [500, 100, 2400, 300]);
plot(t, xP, t, xP_trend);
hold on;
plot(t, yP, t, yP_trend);
title('Säkulare Polbewegung');
xlabel('Zeit (Jahre)');
ylabel('Polkoordinate (arcsec)');
legend('x_P', 'x_P-fit', 'y_P', 'y_P-fit');
saveas(gcf, 'plots/saekulare_polbewegung.png');

% Heatmap
figure;
colormap(hot)
hist3([xP, yP], 'CdataMode', 'auto', 'FaceColor', 'interp');
colorbar;
view(0,90)
title('Heatmap Polkoordinaten');
xlabel('xP (arcsec)');
ylabel('yP (arcsec)');
saveas(gcf, 'plots/heatmap.png');

% Spektrogramm
figure;
[s, f, t_spectro] = spectrogram(detrend(xP), 256, [], [], Fs, 'yaxis');
period = 1 ./ f;
pcolor(t_spectro + t(1), period, abs(s)); % Verschiebung von t_spectro für die korrekte Zeitdarstellung
shading interp;
colorbar;
title('Spektrogramm x_P');
xlabel('Zeit (Jahre)');
ylabel('Periode (Jahre)');
set(gca, 'YScale', 'log'); % Logarithmische Skala für bessere Darstellung
ylim([min(period) max(period)]);
saveas(gcf, 'plots/spektrogramm.png');

% 3D Plot
figure;
set(gcf, 'Position',  [500, 100, 600, 2400]);
plot3(xP, yP, t, '-');
hold on;
plot3([0, 0], [0, 0], [min(t), max(t)], 'r-', 'LineWidth', 2); % Rote Linie bei (0,0)
plot3(polyval(p_xP, t), polyval(p_yP, t), t, 'g-', 'LineWidth', 2); % Säkulare Änderung
grid on;
title('Polbewegung 3D');
xlabel('x_P (arcsec)');
ylabel('y_P (arcsec)');
zlabel('Zeit (Jahre)');
rotate3d on;
view(-40,5);
hold off;

saveas(gcf, 'plots/3d_plot.png');

% Konturplot
[X, Y] = meshgrid(linspace(min(xP), max(xP), 100), linspace(min(yP), max(yP), 80));
Z = griddata(xP, yP, t, X, Y);
figure;
contourf(X, Y, Z);
colorbar;
title('Konturplot Polkoordinaten');
xlabel('xP (arcsec)');
ylabel('yP (arcsec)');
saveas(gcf, 'plots/konturplot.png');

% Histogramm
figure;
subplot(2,1,1);
histogram(xP, 30);
title('Histogramm Polkoordinate x_P');
xlabel('xP (arcsec)');
ylabel('Häufigkeit');

subplot(2,1,2);
histogram(yP, 30);
title('Histogramm der Polkoordinate y_P');
xlabel('yP (arcsec)');
ylabel('Häufigkeit');
saveas(gcf, 'plots/histogramm.png');

% Monatsdurchschnitt über eine Hauptperiode
main_period = min(xP_main_period, yP_main_period); % Hauptperiode bestimmen
months = mod(floor((t - t(1)) * 12 / main_period), 12) + 1; % Monate in der Hauptperiode bestimmen
unique_months = unique(months);
num_months = length(unique_months);
monthly_avg_xP = zeros(num_months, 1);
monthly_avg_yP = zeros(num_months, 1);
for i = 1:num_months
    monthly_avg_xP(i) = mean(xP(months == unique_months(i)));
    monthly_avg_yP(i) = mean(yP(months == unique_months(i)));
end

% Wiederholung über zwei Hauptperioden
repeated_months = [1:num_months num_months+1:2*num_months]'; % Wiederholte Monate
repeated_xP = repmat(monthly_avg_xP, 2, 1);
repeated_yP = repmat(monthly_avg_yP, 2, 1);

figure;
plot(repeated_months, repeated_xP, '-o', 'MarkerFaceColor','auto');
hold on;
plot(repeated_months, repeated_yP, '-o', 'MarkerFaceColor','auto');
title('Monatliche Durchschnittswerte der Polkoordinaten über zwei Hauptperioden');
xlabel('Monat');
ylabel('Polkoordinate (arcsec)');
legend('xP', 'yP');
saveas(gcf, 'plots/monatsdurchschnitt_zwei_perioden.png');

% Continuous Wavelet Transform mit Periodenlänge
figure;
[wt, f] = cwt(xP, 'amor', Fs);
period = 1 ./ f;
pcolor(t, period, abs(wt));
shading interp;
colorbar;
title('Continuous Wavelet Transformation x_P');
xlabel('Zeit (Jahre)');
ylabel('Periode (Jahre)');
set(gca, 'YScale', 'log'); % Logarithmische Skala für bessere Darstellung
ylim([min(period) max(period)]);
saveas(gcf, 'plots/wavelet_transform.png');

% 3D mittel
months = mod(floor((t - t(1)) * 12 / main_period), 12) + 1; % Monate in der Hauptperiode bestimmen
unique_months = unique(months);
num_months = length(unique_months);
monthly_avg_xP = zeros(num_months, 1);
monthly_avg_yP = zeros(num_months, 1);
for i = 1:num_months
    monthly_avg_xP(i) = mean(xP(months == unique_months(i)));
    monthly_avg_yP(i) = mean(yP(months == unique_months(i)));
end
repeated_xP = repmat(monthly_avg_xP, 2, 1);
repeated_yP = repmat(monthly_avg_yP, 2, 1);
repeated_time = (0:(2*num_months-1)) * (main_period / num_months);
figure;
plot3(repeated_time, repeated_xP, repeated_yP, 'o-', 'MarkerFaceColor','auto');
hold on;
plot3(repeated_time, zeros(size(repeated_time)), zeros(size(repeated_time)), 'r-', 'LineWidth', 2); % Rote Linie bei (0,0)
grid on;
title('Polbewegung zwei Hauptperioden');
xlabel('Zeit (Jahre)');
ylabel('xP (arcsec)');
zlabel('yP (arcsec)');
rotate3d on;
hold off;
axis equal
saveas(gcf, 'plots/polbewegung_zwei_perioden.png');


% Ergebnisse ausgeben und in log.txt speichern
logfile = fopen('log.txt', 'w');
fprintf(logfile, 'Hauptperiode in xP: %.4f Jahre\n', xP_main_period);
fprintf(logfile, 'Hauptperiode in yP: %.4f Jahre\n', yP_main_period);
fprintf(logfile, 'Säkulare Änderung in xP: %.4f arcsec/Jahr\n', p_xP(1));
fprintf(logfile, 'Säkulare Änderung in yP: %.4f arcsec/Jahr\n', p_yP(1));
fclose(logfile);

disp(['Hauptperiode in xP: ', num2str(xP_main_period), ' Jahre']);
disp(['Hauptperiode in yP: ', num2str(yP_main_period), ' Jahre']);
disp(['Säkulare Änderung in xP: ', num2str(p_xP(1)), ' arcsec/Jahr']);
disp(['Säkulare Änderung in yP: ', num2str(p_yP(1)), ' arcsec/Jahr']);
