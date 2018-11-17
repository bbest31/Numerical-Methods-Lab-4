%% Part 2 - Brandon Best
load('walk1.mat');
load('human_data.mat');

% M_hip, M_knee, M_foot
M{1,1} = [0.939692000000000,-0.342021000000000,1.06330000000000e-07,0.133632000000000;0.342020000000000,0.939693000000000,-4.99153000000000e-08,-0.175942000000000;-8.28458000000000e-08,8.32722000000000e-08,1,0.0822827000000000;0,0,0,1];
M{1,2} = [1,-7.79704000000000e-17,-2.10284000000000e-08,-1.56828000000000e-06;7.79704000000000e-17,1,1.86594000000000e-15,-0.606984000000000;2.10284000000000e-08,-1.86594000000000e-15,1,5.96046000000000e-07;0,0,0,1];
M{1,3} = [1,6.18174000000000e-08,4.15020000000000e-07,-1.99659000000000e-06;-4.15020000000000e-07,-7.54979000000000e-08,1,-0.553825000000000;6.18174000000000e-08,-1,-7.54979000000000e-08,6.55651000000000e-07;0,0,0,1];
theta = [200,10,10,10]';


[pos,J] = evalRobot3D(M,theta);
% disp(pos);
% disp(J);

for i=1:length(L)
    angleL = invKin3D(Ml,theta,L(:,i));
    angleR = invKin3D(Mr,theta,R(:,i));
    stick3D([angleL;angleR]);
end

function [pos,J] = evalRobot3D(M,theta)
    Rz = [1 0 0 0;...
        0 cos(theta(3)) -sin(theta(3)) 0;...
        0 sin(theta(3)) cos(theta(3)) 0;...
        0 0 0 1];
    Ry = [cos(theta(2)) 0 sin(theta(2)) 0;...
        0 1 0 0;...
        -sin(theta(2)) 0 cos(theta(2)) 0;...
        0 0 0 1];
    
    Rx = [cos(theta(1)) -sin(theta(1)) 0 0;...
        sin(theta(1)) cos(theta(1)) 0 0;...
        0 0 1 0;...
        0 0 0 1];
    RX = [cos(theta(4)) -sin(theta(4)) 0 0;...
        sin(theta(4)) cos(theta(4)) 0 0;...
        0 0 1 0;...
        0 0 0 1];
    pos = M{1,1}*Rz*Ry*Rx*M{1,2}*RX*M{1,3}*[0 0 0 1]';
    P0 = M{1,1}*Rz*[0 0 0 1]';
    P0 = P0(1:3,1);
    P1 = M{1,1}*Rz*Ry*[0 0 0 1]';
    P1 = P1(1:3,1);
    P2 = M{1,1}*Rz*Ry*Rx*M{1,2}*[0 0 0 1]';
    P2 = P2(1:3,1);
    pos = pos(1:3,1);
    J = [P0 P1 P2 pos];

end

function angles = invKin3D(M,theta,pos)
% Newton Method
    [p,J] = evalRobot3D(M,theta);
    x = theta;
    for i=0:10
        f = p - pos;
        s = -J\f;
        x = x + s;
        [p,J] = evalRobot3D(M,x);
    end
    angles = x;
end

function stick3D(theta) 
% Input: theta: angles of each joint, 1x8 vector
% Output: Plot the legs in 3D

M{1,1} = [0.939692000000000,-0.342021000000000,1.06330000000000e-07,0.133632000000000;0.342020000000000,0.939693000000000,-4.99153000000000e-08,-0.175942000000000;-8.28458000000000e-08,8.32722000000000e-08,1,0.0822827000000000;0,0,0,1];
M{1,2} = [1,-7.79704000000000e-17,-2.10284000000000e-08,-1.56828000000000e-06;7.79704000000000e-17,1,1.86594000000000e-15,-0.606984000000000;2.10284000000000e-08,-1.86594000000000e-15,1,5.96046000000000e-07;0,0,0,1];
M{1,3} = [1,6.18174000000000e-08,4.15020000000000e-07,-1.99659000000000e-06;-4.15020000000000e-07,-7.54979000000000e-08,1,-0.553825000000000;6.18174000000000e-08,-1,-7.54979000000000e-08,6.55651000000000e-07;0,0,0,1];
M{1,4} = [0.939693000000000,0.342020000000000,-1.34265000000000e-07,-0.127991000000000;-0.342020000000000,0.939693000000000,-4.49897000000000e-08,-0.175943000000000;1.10781000000000e-07,8.81979000000000e-08,1,0.0822826000000000;0,0,0,1];
M{1,5} = [1,4.42646000000000e-18,5.01038000000000e-09,1.56828000000000e-06;-4.42646000000000e-18,1,4.71169000000000e-16,-0.630394000000000;-5.01038000000000e-09,-4.71169000000000e-16,1,6.55651000000000e-07;0,0,0,1];
M{1,6} = [1,-1.05770000000000e-07,-1.92253000000000e-07,2.11850000000000e-06;1.92253000000000e-07,-7.54979000000000e-08,1,-0.555133000000000;-1.05770000000000e-07,-1,-7.54979000000000e-08,6.55651000000000e-07;0,0,0,1];

% kinematcis
[p_hip_l, p_knee_l, p_foot_l] = kin3D(M(1:3), theta(1:4));
[p_hip_r, p_knee_r, p_foot_r] = kin3D(M(4:6), theta(5:8));

% visualize
vis3D(p_hip_l, p_knee_l, p_foot_l, p_hip_r, p_knee_r, p_foot_r);

end

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

function vis3D(hip_l, knee_l, foot_l, hip_r, knee_r, foot_r)
%% plot the legs in 3D

% get the positions
xl = [hip_l(1) knee_l(1) foot_l(1)];
yl = [hip_l(2) knee_l(2) foot_l(2)];
zl = [hip_l(3) knee_l(3) foot_l(3)];
xr = [hip_r(1) knee_r(1) foot_r(1)];
yr = [hip_r(2) knee_r(2) foot_r(2)];
zr = [hip_r(3) knee_r(3) foot_r(3)];

fig = gcf; 
plot3(xl,yl,zl, 'ro-');
% Set the viewing angle
az = -35;
el = 50;
view(az, el);
axis equal;
axis([-0.5 0.5 -1.5 0 -1 1]);
hold on;
plot3(xr,yr,zr, 'bo-');
% connect two legs
plot3([hip_l(1) hip_r(1)], [hip_l(2) hip_r(2)], [hip_l(3) hip_r(3)], 'k-');

xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
hold off;
legend('Left Leg','Right Leg','Hip');
drawnow();
end

