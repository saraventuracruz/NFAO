% Define the input grid
[x, y] = meshgrid(-10: 10);
%% Calculate the two surfaces
z1 = 3-x;
z2 = -y;
%% Visualize the two surfaces
surf(x, y, z1, 'FaceColor', [0.8 0.8 0.8]); 
hold on
surf(x, y, z2, 'FaceColor', [0.3 0.3 0.3]);
view(3);
%% Take the difference between the two surface heights and find the contour
% where that surface is zero.
zdiff = z1 - z2;
C = contours(x, y, zdiff, [0 0]);
%% Extract the x- and y-locations from the contour matrix C.
xL = C(1, 2:end);
yL = C(2, 2:end);
% Interpolate on the first surface to find z-locations for the intersection
% line.
zL = interp2(x, y, z1, xL, yL);
% Visualize the line.
line(xL, yL, zL, 'Color', 'k', 'LineWidth', 3);