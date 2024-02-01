function [] = plotEarthOrbit(r_V, R, figTitle)
%PLOTORBIT Summary of this function goes here
%   Detailed explanation goes here

% get coordinates of a sphere that represents earth
[surfX, surfY, surfZ] = sphere;
surfX = surfX*R;
surfY = surfY*R;
surfZ = surfZ*R;

figure
plot3(r_V(1,:), r_V(2,:), r_V(3,:), 'LineWidth', 1.5)
hold on
surf(surfX, surfY, surfZ)
axis equal
title(figTitle)
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on

end

