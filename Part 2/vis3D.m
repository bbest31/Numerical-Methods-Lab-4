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
