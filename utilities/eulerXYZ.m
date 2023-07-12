function R = eulerXYZ(a1, a2, a3)
% Convert XYZ Euler angles to rotation matrix

R1 = [
1, 0, 0;
0, cos(a1), -sin(a1);
0, sin(a1), cos(a1)];

R2 = [
cos(a2), 0, sin(a2);
0, 1, 0;
-sin(a2), 0, cos(a2)];

R3 = [
cos(a3), -sin(a3), 0;
sin(a3), cos(a3), 0;
0, 0, 1];

R = R1*R2*R3;

end