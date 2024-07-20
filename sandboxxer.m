clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLES

% Initialize the simulation with Neumayer Station's coordinates
phi = -70.6734; % latitude in degrees
lambda = -8.2741; % longitude in degrees
% % phi = 85; % latitude in degrees
% % lambda = 90; % longitude in degrees
% phi = 48; % latitude in degrees
% lambda = -7; % longitude in degrees

% Date
yr = 2024;
mt = 9;
dy = 21;
ut1 = 12;

[R0,P0, N0, G0, R1,P1, N1, G1,t, mjd, jd, gast, gmst, eps0, deleps, delpsi, za, thetaa, zetaa] = GDS_INT_TO_LCL(lambda, phi, yr, mt, dy, ut1);

figure;
colormap('gray');
set(gcf, 'Position',  [1720, -7, 1720, 1440]);
tiledlayout(3,3);
nexttile;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [CELESTIAL] CONVENTIONAL INERTIAL [NCP]
hold on; % We want to plot all at once
view(30,30); % Good perspective setting
axis equal; % Isoview
title('Inertial System i_0')

% [1] Make a sphere, plot it
[u, v, w] = ellipsoid(0,0,0,1,1,1,20);
s = surf(u, v, w, 'EdgeColor', 'interp', 'FaceAlpha', 0.2);

% [2] Create the GC and plot them
theta = linspace(-pi,pi,100); % Generate the angles
x = cos(theta); 
y = sin(theta); 
z = zeros(size(theta)); % Process to circles

gc_xy_i = [x' y' z']';
gc_xz_i = [x' z' y']';
gc_yz_i = [z' y' x']';

plot3(gc_xy_i(1,:), gc_xy_i(2,:), gc_xy_i(3,:), 'k', 'LineWidth', 1);
plot3(gc_xz_i(1,:), gc_xz_i(2,:), gc_xz_i(3,:), 'k', 'LineWidth', 1);
plot3(gc_yz_i(1,:), gc_yz_i(2,:), gc_yz_i(3,:), 'k', 'LineWidth', 1);

% [5] Create the Ecliptic (black)
gc_ec_i = rotx(eps0)*N1*P1*gc_xy_i; % Positive, cause we rotate the object!
plot3(gc_ec_i(1,:), gc_ec_i(2,:), gc_ec_i(3,:), 'k', 'LineWidth', 2);

% [3] Create the GWC (red semicircle) and plot it
theta = theta/2;
x = cos(theta); y = zeros(size(theta)); z = sin(theta);
gw_i = [x' y' z']';
plot3(gw_i(1,:), gw_i(2,:), gw_i(3,:), 'r', 'LineWidth', 2);

% [6] Create the Ecliptic and plot it
gw_c = G1*N1*P1*gw_i;
h = plot3(gw_c(1,:), gw_c(2,:), gw_c(3,:), 'b', 'LineWidth', 2);




% view(0,90)

% [4] Create inertial transformation matrix and rotate GWC
% for j = 1:0.1:24
%     [R0,P0, N0, G0, R1,P1, N1, G1,t, mjd, jd, gast, gmst, eps0, deleps, delpsi, za, thetaa, zetaa] = GDS_INT_TO_LCL(lambda, phi, yr, mt, dy, j);
%     gw_c = G1*N1*P1*gw_i;
%     set(h, 'XData', gw_c(1, :), 'YData', gw_c(2, :), 'ZData', gw_c(3, :));
%     if j == 1
%         pause(1)
%     end
%     pause(0.05);
% end

[x y z] = sph2cart(deg2rad(lambda),deg2rad(phi),1);
POS0 = [x y z]';
POS1 = G1*N1*P1*POS0;

K0 = 0.5*eye(3);
K1 = R0*K0;
quiver3(0,0,0, K0(1,1), K0(1,2), K0(1,3), 'r', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);
quiver3(0,0,0, K0(2,1), K0(2,2), K0(2,3), 'g', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);
quiver3(0,0,0, K0(3,1), K0(3,2), K0(3,3), 'b', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);

K0 = K1;
quiver3(POS1(1),POS1(2),POS1(3), K0(1,1), K0(1,2), K0(1,3), 'r', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);
quiver3(POS1(1),POS1(2),POS1(3), K0(2,1), K0(2,2), K0(2,3), 'g', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);
quiver3(POS1(1),POS1(2),POS1(3), K0(3,1), K0(3,2), K0(3,3), 'b', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);

quiver3(0,0,0, POS1(1),POS1(2),POS1(3), 'm', 'filled','LineWidth',2, 'AutoScale','off','MaxHeadSize',0.5);
% quiver3(0,0,0, 0, 0, 0, 'g', 'filled','LineWidth',2);
text(0,0,1.3,'Local Time UT1: '+ string(ut1)+'h')
% plot3()
% [4] Create the Sunposition and plot it
[aapp,dapp,~] = GDS_SOLARPOS(t);
[x y z] = sph2cart(deg2rad(aapp),deg2rad(dapp),1);

quiver3(POS1(1),POS1(2),POS1(3),x,y,z, 'y', 'filled','LineWidth',2, 'AutoScale','off','MaxHeadSize',0.5);

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [CELESTIAL] CONVENTIONAL INERTIAL
nexttile
hold on; % We want to plot all at once
view(30,30); % Good perspective setting
axis equal; % Isoview
title('Inverse Test')

% [1] Make a sphere, plot it
[u, v, w] = ellipsoid(0,0,0,1,1,1,20);
s = surf(u, v, w, 'EdgeColor', 'interp', 'FaceAlpha', 0.2);


M = inv(R0);
K0 = M*K0;
POSS = R0*[x y z]';
quiver3(0,0,0, K0(1,1), K0(1,2), K0(1,3), 'r', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);
quiver3(0,0,0, K0(2,1), K0(2,2), K0(2,3), 'g', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);
quiver3(0,0,0, K0(3,1), K0(3,2), K0(3,3), 'b', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);
quiver3(0,0,0,POSS(1),POSS(2),POSS(3),'y', 'filled','LineWidth',2, 'AutoScale','off','MaxHeadSize',0.5);

hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [CELESTIAL] CONVENTIONAL INERTIAL
nexttile
hold on; % We want to plot all at once
view(30,30); % Good perspective setting
axis equal; % Isoview
title('Mean inertial')

% [1] Make a sphere, plot it
[u, v, w] = ellipsoid(0,0,0,1,1,1,20);
s = surf(u, v, w, 'EdgeColor', 'flat', 'FaceAlpha', 0.3);


gw_c = G1*gw_i;
h = plot3(gw_c(1,:), gw_c(2,:), gw_c(3,:), 'b', 'LineWidth', 2);

[x y z] = sph2cart(deg2rad(lambda),deg2rad(phi),1);
POS0 = [x y z]';
POS1 = G1*POS0;

% Originales Koordinatensystem
K0 = 0.5*eye(3);
quiver3(0,0,0, K0(1,1), K0(1,2), K0(1,3), 'r', 'filled','LineWidth',3, 'AutoScale','off');
quiver3(0,0,0, K0(2,1), K0(2,2), K0(2,3), 'g', 'filled','LineWidth',3, 'AutoScale','off');
quiver3(0,0,0, K0(3,1), K0(3,2), K0(3,3), 'b', 'filled','LineWidth',3, 'AutoScale','off');

% Nach kompletter Rotation und Translation 
K1 = R0*K0;
quiver3(POS1(1),POS1(2),POS1(3), K1(1,1), K1(1,2), K1(1,3), 'r', 'filled','LineWidth',3, 'AutoScale','off');
quiver3(POS1(1),POS1(2),POS1(3), K1(2,1), K1(2,2), K1(2,3), 'g', 'filled','LineWidth',3, 'AutoScale','off');
quiver3(POS1(1),POS1(2),POS1(3), K1(3,1), K1(3,2), K1(3,3), 'b', 'filled','LineWidth',3, 'AutoScale','off');

quiver3(0,0,0, POS1(1),POS1(2),POS1(3), 'm', 'filled','LineWidth',2, 'AutoScale','off');
% quiver3(0,0,0, 0, 0, 0, 'g', 'filled','LineWidth',2);

% [1] Make a sphere, plot it
[u, v, w] = ellipsoid(POS1(1),POS1(2),POS1(3),.3,.3,.3,10);
s = surf(u, v, w, 'EdgeColor', 'flat', 'FaceAlpha', 0.3);
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
% sunpos(:,:,2) = zeros(1,2);

for mt = 1:1:6
    % ut1 = 0;
    % [aapp, dapp, ~] = GDS_SOLARPOS(t);
    % [x, y, z] = sph2cart(deg2rad(aapp), deg2rad(dapp), 1);
    % vec = R0 * [x y z]';
    % [az, el, ~] = cart2sph(vec(1), vec(2), vec(3));
    % sunpos = [rad2deg(az) rad2deg(el)];
    i = 1;
    for ut1 = 0.0:0.1:24
        % Acquiring exact position per hour in local coordinates
        [R0,P0, N0, G0, R1,P1, N1, G1,t, mjd, jd, gast, gmst, eps0, deleps, delpsi, za, thetaa, zetaa] = GDS_INT_TO_LCL(lambda, phi, yr, mt, dy, ut1);
        [aapp, dapp, ~] = GDS_SOLARPOS(t); % Position of the sun in inertial coordinates
        [x, y, z] = sph2cart(deg2rad(aapp), deg2rad(dapp), 1); % Conversion to cartesian coordinates
        vec = R0 * [x y z]';  % Transformation to local coordinates
        [az, el, ~] = cart2sph(vec(1), vec(2), vec(3)); % Conversion back to spherical coordinates
        
        sunpos(i, :) = [rad2deg(az) rad2deg(el) ut1];
        i = i+1;
    end
    pos_data(:,:,dy,mt) = sunpos;

    sunpos = sunpos(sunpos(:, 2) >= 0, :); % Filter out negative elevations
    % Plot and store handle
    hsky = GDS_SKYPLOT(sunpos(:, 1), sunpos(:, 2), '-', 3, colors(mt, :), '');
    if length(sunpos)>1
        legend_handles(end+1) = hsky; % Collect handles
        legend_entries{end+1} = month_names{mt}; % Collect legend entries
    end
    % Transform data to Cartesian coordinates
    yy = (90 - sunpos(sunpos(:,3)==12, 2)) .* cos(sunpos(sunpos(:,3)==12, 1) / 180 * pi);
    xx = (90 - sunpos(sunpos(:,3)==12, 2)) .* sin(sunpos(sunpos(:,3)==12, 1) / 180 * pi);
    plot(xx,yy,'MarkerSize',9,'Color','r','Marker','o','MarkerEdgeColor','k', 'MarkerFaceColor','y');
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

for mt = 7:1:12
    % ut1 = 0;
    % [aapp, dapp, ~] = GDS_SOLARPOS(t);
    % [x, y, z] = sph2cart(deg2rad(aapp), deg2rad(dapp), 1);
    % vec = R0 * [x y z]';
    % [az, el, ~] = cart2sph(vec(1), vec(2), vec(3));
    % sunpos = [rad2deg(az) rad2deg(el)];
    i = 1;
    for ut1 = 0.0:0.1:24
        % Acquiring exact position per hour in local coordinates
        [R0,P0, N0, G0, R1,P1, N1, G1,t, mjd, jd, gast, gmst, eps0, deleps, delpsi, za, thetaa, zetaa] = GDS_INT_TO_LCL(lambda, phi, yr, mt, dy, ut1);
        [aapp, dapp, ~] = GDS_SOLARPOS(t); % Position of the sun in inertial coordinates
        [x, y, z] = sph2cart(deg2rad(aapp), deg2rad(dapp), 1); % Conversion to cartesian coordinates
        vec = R0 * [x y z]';  % Transformation to local coordinates
        [az, el, ~] = cart2sph(vec(1), vec(2), vec(3)); % Conversion back to spherical coordinates
        
        sunpos(i, :) = [rad2deg(az) rad2deg(el) ut1];
        i = i+1;
    end
    pos_data(:,:,dy,mt) = sunpos;

    sunpos = sunpos(sunpos(:, 2) >= 0, :); % Filter out negative elevations
    % Plot and store handle
    hsky = GDS_SKYPLOT(sunpos(:, 1), sunpos(:, 2), '-', 3, colors(mt, :), '');
    % set(hsky,'XData',0, 'YData',0)
    if length(sunpos)>1
        legend_handles(end+1) = hsky; % Collect handles
        legend_entries{end+1} = month_names{mt}; % Collect legend entries
    end
    % Transform data to Cartesian coordinates
    yy = (90 - sunpos(sunpos(:,3)==12, 2)) .* cos(sunpos(sunpos(:,3)==12, 1) / 180 * pi);
    xx = (90 - sunpos(sunpos(:,3)==12, 2)) .* sin(sunpos(sunpos(:,3)==12, 1) / 180 * pi);
    plot(xx,yy,'MarkerSize',9,'Color','r','Marker','o','MarkerEdgeColor','k', 'MarkerFaceColor','y')
end

% Update legend at the end
legend(legend_handles, legend_entries);
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL
nexttile
nexttile;
hold on;
view(110,10);
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
% theta = theta/2;
% x = cos(theta); y = zeros(size(theta)); z = sin(theta);
% gw_i = [x' y' z']';
% plot3(gw_i(1,:), gw_i(2,:), gw_i(3,:), 'r', 'LineWidth', 2);

% [3] Create the Spring Point and plot it
% [t, ~, ~] = GDS_JULIANC(yr, mt, dy, ut1);
% [aapp,dapp,~] = GDS_SOLARPOS(t);
% [x, y, z] = sph2cart(deg2rad(aapp),deg2rad(dapp),1);
% pos = POINT_TO_VECTOR(x,y,z);
% HOLD_PLOT_VEC3(pos, 'Spring Point: '+string(ut1)+'h');

% [4] Create inertial transformation matrix and rotate GWC
% [R0, G0] = GDS_INT_TO_LCL(lambda,phi,yr,mt,dy,ut1);
% gw_i = G0*gw_i;
% plot3(gw_i(1,:), gw_i(2,:), gw_i(3,:), 'b', 'LineWidth', 2);

% [5] Create the Ecliptic and plot it
% gc_ec_i = rotx(23.5)*gc_xy_i;
% plot3(gc_ec_i(1,:), gc_ec_i(2,:), gc_ec_i(3,:), 'k', 'LineWidth', 2);

ut1 = 12;
K0 = 0.5*eye(3);

quiver3(0,0,0, K0(1,1), K0(1,2), K0(1,3), 'r', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);
quiver3(0,0,0, K0(2,1), K0(2,2), K0(2,3), 'g', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);
quiver3(0,0,0, K0(3,1), K0(3,2), K0(3,3), 'b', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);
for k = 1:1:6
    sunpos = pos_data(:,:,dy,k);
    sunpos(:,1:2) = deg2rad(sunpos(:,1:2));
    r = ones(size(sunpos(:,1)));
    [x y z] = sph2cart(sunpos(:, 1), sunpos(:, 2), r);
    plot3(x,y,z,'Color',colors(k, :), 'LineWidth', 2);


    [x y z] = sph2cart(sunpos(sunpos(:,3)==12, 1), sunpos(sunpos(:,3)==12, 2), 1);
    ut1 = 12;
    if z
        % quiver3(0,0,0,x,y,z,'y', 'filled','LineWidth',2, 'AutoScale','off','MaxHeadSize',0.5);
        plot3(x,y,z,'MarkerSize',10,'Color','y','Marker','o','MarkerEdgeColor','k', 'MarkerFaceColor','y');

    end
    mday = pos_data(pos_data(:, 2, dy, k) >= 0, :, dy, k);
    if mday
        disp('Sonnenaufgang für '+string(dy)+'/'+month_names(k)+'/'+string(yr) + ': '+string(floor(mday(1,3)))+ 'h ' + string((mday(1,3)-floor(mday(1,3)))*60)+'min');
        disp('Sonnenuntergang für  '+string(dy)+'/'+month_names(k)+'/'+string(yr) + ': '+string(floor(mday(end,3)))+ 'h ' + string((mday(end,3)-floor(mday(end,3)))*60)+'min');
    else
        disp('Kein Sonnenauf/ -untergang für: '+string(dy)+'/'+month_names(k)+'/'+string(yr))
    end
end
zlim([0 1]);

hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL

nexttile;
hold on;
view(110,10);
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
% theta = theta/2;
% x = cos(theta); y = zeros(size(theta)); z = sin(theta);
% gw_i = [x' y' z']';
% plot3(gw_i(1,:), gw_i(2,:), gw_i(3,:), 'r', 'LineWidth', 2);

% [3] Create the Spring Point and plot it
% [t, ~, ~] = GDS_JULIANC(yr, mt, dy, ut1);
% [aapp,dapp,~] = GDS_SOLARPOS(t);
% [x, y, z] = sph2cart(deg2rad(aapp),deg2rad(dapp),1);
% pos = POINT_TO_VECTOR(x,y,z);
% HOLD_PLOT_VEC3(pos, 'Spring Point: '+string(ut1)+'h');

% [4] Create inertial transformation matrix and rotate GWC
% [R0, G0] = GDS_INT_TO_LCL(lambda,phi,yr,mt,dy,ut1);
% gw_i = G0*gw_i;
% plot3(gw_i(1,:), gw_i(2,:), gw_i(3,:), 'b', 'LineWidth', 2);

% [5] Create the Ecliptic and plot it
% gc_ec_i = rotx(23.5)*gc_xy_i;
% plot3(gc_ec_i(1,:), gc_ec_i(2,:), gc_ec_i(3,:), 'k', 'LineWidth', 2);
ut1 = 12;
K0 = 0.5*eye(3);

quiver3(0,0,0, K0(1,1), K0(1,2), K0(1,3), 'r', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);
quiver3(0,0,0, K0(2,1), K0(2,2), K0(2,3), 'g', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);
quiver3(0,0,0, K0(3,1), K0(3,2), K0(3,3), 'b', 'filled','LineWidth',3, 'AutoScale','off','MaxHeadSize',0.5);
for k = 7:1:12
    sunpos = pos_data(:,:,dy,k);
    sunpos(:,1:2) = deg2rad(sunpos(:,1:2));
    r = ones(size(sunpos(:,1)));
    [x y z] = sph2cart(sunpos(:, 1), sunpos(:, 2), r);
    plot3(x,y,z,'Color',colors(k, :), 'LineWidth', 2);


    [x y z] = sph2cart(sunpos(sunpos(:,3)==12, 1), sunpos(sunpos(:,3)==12, 2), 1);
    ut1 = 12;
    if z
        % quiver3(0,0,0,x,y,z,'y', 'filled','LineWidth',2, 'AutoScale','off','MaxHeadSize',0.5);
        plot3(x,y,z,'MarkerSize',10,'Color','y','Marker','o','MarkerEdgeColor','k', 'MarkerFaceColor','y');
    end

    mday = pos_data(pos_data(:, 2, dy, k) >= 0, :, dy, k);
    if mday
        disp('Sonnenaufgang für '+string(dy)+'/'+month_names(k)+'/'+string(yr) + ': '+string(floor(mday(1,3)))+ 'h ' + string((mday(1,3)-floor(mday(1,3)))*60)+'min');
        disp('Sonnenuntergang für  '+string(dy)+'/'+month_names(k)+'/'+string(yr) + ': '+string(floor(mday(end,3)))+ 'h ' + string((mday(end,3)-floor(mday(end,3)))*60)+'min');
    else
        disp('Kein Sonnenauf/ -untergang für: '+string(dy)+'/'+month_names(k)+'/'+string(yr))
    end
end




zlim([0 1]);

hold off;
nexttile