function [] = plotBody(position, radius, color)

% get coordinates of a sphere that represents earth
[surfX, surfY, surfZ] = sphere;
surfX = surfX*radius;
surfY = surfY*radius;
surfZ = surfZ*radius;
x = position(1);
y = position(2);
z = position(3);
surf(surfX+x, surfY+y, surfZ+z, 'FaceColor', color, 'EdgeColor','none')

end

