function [vertices, sideFaces, bottomFaces] = calcCylinder(r, R, Radius, Height, n)
% Return vertices and faces of a Cylinder geometry
% n SideCount
% R Rotation from parent coordinate frame
% r Translation from parent coordinate frame

% Vertices
vertices_0 = zeros(2*n, 3);
for i = 1:n
    theta = 2*pi/n*(i-1);
    vertices_0(i,:) = [Radius*cos(theta), Radius*sin(theta), -Height/2];
    vertices_0(n+i,:) = [Radius*cos(theta), Radius*sin(theta), Height/2];
end

vertices = r' + vertices_0*R';

% Side faces
sideFaces = zeros(n, 4);
for i = 1:(n-1)
    sideFaces(i,:) = [i, i+1, n+i+1, n+i];
end
sideFaces(n,:) = [n, 1, n+1, 2*n];

% Bottom faces
bottomFaces = [
    1:n;
    (n+1):2*n];

end