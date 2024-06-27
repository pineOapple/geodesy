function INSERT_GREAT_CRICLE(AX, AY, AZ)
    % Define the circle in the xy-plane
    theta = linspace(0, 2*pi, 100);
    x = cos(theta);
    y = sin(theta);
    z = zeros(size(theta));
    
    % Create rotation matrices
    Rx = rotx(deg2rad(AX));
    Ry = roty(deg2rad(AY));
    Rz = rotz(deg2rad(AZ));
    
    % Combine rotations (adjust sequence as necessary for specific rotation order)
    R = Rz * Ry * Rx;  % For example, rotate about z, then y, then x
    
    % Apply rotation to circle coordinates
    % Points are in rows: transpose them for matrix multiplication
    points = [x; y; z];
    rotatedPoints = R * points;
    
    % Extract rotated coordinates
    xRot = rotatedPoints(1, :);
    yRot = rotatedPoints(2, :);
    zRot = rotatedPoints(3, :);
    
    % Plot the original and rotated circles
    plot3(x, y, z, 'k', 'LineWidth', 1);  % Original circle in xy-plane
    plot3(x, z, y, 'k', 'LineWidth', 1);  % Original circle in xz-plane
    plot3(z, y, x, 'k', 'LineWidth', 1);  % Original circle in yz-plane
    
    plot3(xRot, yRot, zRot, 'r', 'LineWidth', 1);  % Rotated circle
    
    % Set labels and title
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(sprintf('Circle rotated by X: %d°, Y: %d°, Z: %d°', angleX, angleY, angleZ));
    
    hold off;
end
