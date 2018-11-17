%% Lab 4 - Part 1 - Brandon Best

links = [1.5;1];
theta = [45;45];
a = 0.005;

[pos,J] = evalRobot2D(links,theta);
% disp(pos);
% disp(J);

J2 = fdJacob2D(links,theta,a);
% disp(J2);

% plotRobot2D(links,theta);

testInvKin(10);

function testInvKin(n)
    ls = [0.5,0.5]';
    t = rand(2,1);
    clf;
    plotRobot2D(ls,t);
    hold off;
    while(1)
        desired=ginput(1)';
    
        clf;
        plot(desired(1),desired(2),'*');
        hold on;
        plotRobot2D(ls,t,':');
    
        t = invKin2D(ls,t,desired,n,0);
        plotRobot2D(ls,t);
        hold off;
    end
end

% Function takes in link lengths and theta and returns the end effector
% position and Jacobian after rotating the joints by theta.
function [pos,J] = evalRobot2D (l,theta)
    link_1 = l(1,1);
    link_2 = l(2,1);
    theta_1 = theta(1,1);
    theta_2 = theta(2,1);
    
    link_1_pos = [link_1*cos(theta_1) link_1*sin(theta_1)];
    link_2_pos = [link_1*cos(theta_1)+link_2*cos(theta_1+theta_2) link_1*sin(theta_1)+link_2*sin(theta_1+theta_2)];
    pos = link_2_pos';
    
    f1 = [-link_1*sin(theta_1)-link_2*sin(theta_1+theta_2) -link_2*sin(theta_1+theta_2)];
    f2 = [link_1*cos(theta_1)+link_2*cos(theta_1+theta_2) link_2*cos(theta_1+theta_2)];
    J = [f1;f2];
    return;
end



function J = fdJacob2D(l, theta, alpha)

J = (evalRobot2D(l,theta+[alpha;0])-evalRobot2D(l,theta-[alpha;0]))/(2*alpha);
J = [J',((evalRobot2D(l,theta+[0;alpha])-evalRobot2D(l,theta-[0;alpha]))/(2*alpha))'];
end

function angles = invKin2D(l,theta,pos,n,mode)

if(mode == 1)
    % Broyden Update
    [p,J] = evalRobot2D(l,theta);
    delta_q = [];
    e = [];
    for i=0:n
        e = pos - p;
        delta_q = inv(J)*e;
        numerator = p-J*delta_q;
        denominator = delta_q'*delta_q;
        J = J + (numerator/denominator)*delta_q';
        [p,j] = evalRobot2D(l,delta_q);
    end
    
    angles = delta_q;
    
elseif(mode == 0)
    % Newton Method
    [p,J] = evalRobot2D(l,theta);
    x = theta;
    for i=0:n
        f = p - pos;
        s = -J\f;
        x = x + s;
        [p,J] = evalRobot2D(l,x);
    end
    angles = x;
end

end


% plotRobot2D(ls,angs)
% ls: a vector containing the link lengths
% angs: a vector containing the link rotation angles
%
function plotRobot2D(ls,angs,linespec)
if(nargin<3)
  linespec='';
end
holdon=ishold;
hold on;
axis equal;

%draw invisible workspace
ltot=sum(ls);
plot([-ltot,ltot,ltot,-ltot],[-ltot,-ltot,ltot,ltot],'w');

%draw each segment
f=0;
t=0;
for j=1:length(ls)
  t=t+angs(j);
  p=f+(cos(t)+i*sin(t))*(ls(j));
  plotSeg(f,p,linespec);
  hold on;
  f=p;
end

%Restore hold
if holdon
  hold on
else
  hold off
end
end
function plotSeg(from,to,linespec)
d=to-from;
l=norm(d);
r1=l/10;
r2=l/15;
r3=l/7;
d=d/l;
ks=strcat('k',linespec);

mycircle(l/10,real(from),imag(from),ks);
mycircle(l/15,real(to),imag(to),ks);

pts=[from+r1*i*d,to+r2*i*d];
plot(real(pts),imag(pts),ks);
pts=[from-r1*i*d,to-r2*i*d];
plot(real(pts),imag(pts),ks);

%The axis
pts=[to,to+r3*d];
plot(real(pts),imag(pts),strcat('r',linespec));
pts=[to,to+i*r3*d];
plot(real(pts),imag(pts),strcat('g',linespec));
end

function mycircle(rad,px,py,opts)
t=linspace(0,2*pi,10);
plot(rad*cos(t)+px,rad*sin(t)+py,opts);
end
