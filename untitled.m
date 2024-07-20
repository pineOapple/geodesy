% Step 1: Generate the ellipsoid coordinates
[u, v, w] = ellipsoid(0, 0, 0, 1, 1, 1, 20);

% Step 2: Flatten the coordinates into a list of points
numPoints = numel(u);
points = [u(:)'; v(:)'; w(:)'];

% Step 3: Define a custom rotation matrix (3x3)
theta = pi; % Example rotation angle (30 degrees)
R = [cos(theta) 0 -sin(theta); 
    0 1 0;
     sin(theta) 0 cos(theta)]; % Rotation matrix around the z-axis

% Step 4: Apply the rotation matrix to each point
rotatedPoints = R * points;

% Step 5: Reshape the rotated points back to the ellipsoid coordinate shape
u_rot = reshape(rotatedPoints(1, :), size(u));
v_rot = reshape(rotatedPoints(2, :), size(v));
w_rot = reshape(rotatedPoints(3, :), size(w));

% Plot the original ellipsoid
figure;
subplot(1, 2, 1);
surf(u, v, w, 'EdgeColor', 'none');
axis equal;
title('Original Ellipsoid');
xlabel('X');
ylabel('Y');
zlabel('Z');

% Plot the rotated ellipsoid
subplot(1, 2, 2);
surf(u_rot, v_rot, w_rot, 'EdgeColor', 'none');
axis equal;
title('Rotated Ellipsoid');
xlabel('X');
ylabel('Y');
zlabel('Z');
