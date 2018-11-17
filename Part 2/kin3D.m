function [p_hip, p_knee, p_foot] = kin3D(M, theta)
% returns the positions of the links
% theta is a 4x1 vector

Rz3 = [cos(theta(3)) -sin(theta(3)) 0 0; -sin(theta(3)) cos(theta(3)) 0 0; 0 0 1 0; 0 0 0 1];
Ry2 = [cos(theta(2)) 0 sin(theta(2)) 0; 0 1 0 0;-sin(theta(2)) 0 cos(theta(2)) 0; 0 0 0 1];
Rx1 = [1 0 0 0; 0 cos(theta(1)) -sin(theta(1)) 0;0 sin(theta(1)) cos(theta(1)) 0; 0 0 0 1];
Rx4 = [1 0 0 0; 0 cos(theta(4)) -sin(theta(4)) 0;0 sin(theta(4)) cos(theta(4)) 0; 0 0 0 1];

hip = M{1} * [0,0,0,1]';
p_hip = hip(1:3);

knee = M{1} * Rz3 * Ry2 * Rx1 * M{2} * [0,0,0,1]';
p_knee = knee(1:3);

foot = M{1} * Rz3 * Ry2 * Rx1 * M{2} * Rx4 * M{3} * [0,0,0,1]';
p_foot = foot(1:3);
end