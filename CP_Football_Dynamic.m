% This Program is to simulate the realistic football trajectory numerically 
% using Euler method. The accuracy of this program is in one order. The user
% just need to assign some initial conditions to start the program.

clc
clear


% Assigning the constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 10000;
j = 0;
dt = 0.001;

R = 0.11;
A = pi*R^2;

rho = 1.225;
C_D = 0.30; %0.30
C_L = input('Enter a number for lift coefficient (0 - 0.26): ');
C_DM = 0.05;

m = 0.4;
J = (2/5)*m*R^2;
g = -9.81;

k = (1/2)*rho*A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Assigning functions for ODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_wx = @(wx,wy,wz) -(C_DM/J)*k*wx/sqrt(wx^2+wy^2+wz^2);
f_wy = @(wx,wy,wz) -(C_DM/J)*k*wy/sqrt(wx^2+wy^2+wz^2);
f_wz = @(wx,wy,wz) -(C_DM/J)*k*wz/sqrt(wx^2+wy^2+wz^2);

f_x = @(vx,vy,vz,wx,wy,wz) -(C_D/m)*k*vx*sqrt(vx^2+vy^2+vz^2)+(C_L/m)*k*(vx^2+vy^2+vz^2)*(wy*vz-wz*vy)*(1/norm(cross([wx wy wz],[vx vy vz])));
f_y = @(vx,vy,vz,wx,wy,wz) -(C_D/m)*k*vy*sqrt(vx^2+vy^2+vz^2)+(C_L/m)*k*(vx^2+vy^2+vz^2)*(wz*vx-wx*vz)*(1/norm(cross([wx wy wz],[vx vy vz])));
f_z = @(vx,vy,vz,wx,wy,wz) -(C_D/m)*k*vz*sqrt(vx^2+vy^2+vz^2)+(C_L/m)*k*(vx^2+vy^2+vz^2)*(wx*vy-wy*vx)*(1/norm(cross([wx wy wz],[vx vy vz])))+g;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate the arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = zeros(n,1);

x = zeros(n,1);
y = zeros(n,1);
z = zeros(n,1);

wx = zeros(n,1);
wy = zeros(n,1);
wz = zeros(n,1);

vx = zeros(n,1);
vy = zeros(n,1);
vz = zeros(n,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prompt user to enter the informations about the system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Starting point of the ball. Type the following for specific position. \n \n');
fprintf('1 for corner \n');
fprintf('2 for common free kick position \n');
fprintf('3 for center \n');
fprintf('other letter for self define initial position \n \n');
init_choice = input('Please type: ');

    % User choice on the position of football.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch init_choice
    case 1
        x(1) = -30;
        y(1) = 50;
        z(1) = 0;
    case 2
        x(1) = -22;
        y(1) = 30;
        z(1) = 0;
    case 3
        x(1) = 0;
        y(1) = 0;
        z(1) = 0;
    otherwise
        x(1) = input('initial x:');
        y(1) = input('initial y:');
        z(1) = input('initial z:');
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = (pi/180)*input('\nEnter azimuth angle (-180 < a < 180): ');
b = (pi/180)*input('\nEnter elevation angle (0 < b < 180) (usually around 25): ');


v_init = input('\nInitial velocity of ball (usually around 20 - 30 m/s): ');

fprintf('\nInput the angular velocity, wx, wy, wz. One of them must be non zero. (Can be 0.000001)  \n');
wx(1) = input('Initial angular velocity of ball, wx: ');
wy(1) = input('Initial angular velocity of ball, wy: ');
wz(1) = input('Initial angular velocity of ball, wz: ');

vp = input('Enter the view point: ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Assigning initial condition for velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vx(1) = v_init*cos(b)*sin(a);
vy(1) = v_init*cos(b)*cos(a);
vz(1) = v_init*sin(b);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Looping for RK4(5) algorithms n times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n
    wx(i+1) = wx(i)+dt*f_wx(wx(i),wy(i),wz(i));
    wy(i+1) = wy(i)+dt*f_wy(wx(i),wy(i),wz(i));
    wz(i+1) = wz(i)+dt*f_wz(wx(i),wy(i),wz(i));

    vx(i+1) = vx(i)+dt*f_x(vx(i),vy(i),vz(i),wx(i),wy(i),wz(i));
    vy(i+1) = vy(i)+dt*f_y(vx(i),vy(i),vz(i),wx(i),wy(i),wz(i));
    vz(i+1) = vz(i)+dt*f_z(vx(i),vy(i),vz(i),wx(i),wy(i),wz(i));

    x(i+1) = x(i)+vx(i)*dt;
    y(i+1) = y(i)+vy(i)*dt;
    z(i+1) = z(i)+vz(i)*dt;

    t(i+1) = t(i)+dt;

    if x(i+1) < -35 || x(i+1) > 35 || z(i+1) < 0 || y(i+1) < -55 || y(i+1) > 55
        j = i;
        break
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Delete all the unnecessary information such as the position outside the 
% football field from arrays.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x(j+1:end) = [];
y(j+1:end) = [];
z(j+1:end) = [];
wx(j+1:end) = [];
wy(j+1:end) = [];
wz(j+1:end) = [];
vx(j+1:end) = [];
vy(j+1:end) = [];
vz(j+1:end) = [];
t(j+1:end) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Assign coordinate for drawing football field.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = sqrt(vx.^2+vy.^2+vz.^2);
w = sqrt(wx.^2+wy.^2+wz.^2);

Rangex = -9.15:0.1:9.15;
Cy1 = sqrt(9.15^2-Rangex.^2);
Cy2 = -sqrt(9.15^2-Rangex.^2);

Range2x = -7.45:0.1:7.45;
Sy1 = sqrt(9.15^2-Range2x.^2)-39;
Sy2 = -sqrt(9.15^2-Range2x.^2)+39;

corx = 29:0.1:30;
cory = sqrt(1-(corx-30).^2)-50;

cor2x = -30:0.1:-29;
cor2y = sqrt(1-(cor2x+30).^2)-50;

cor3x = -30:0.1:-29;
cor3y = -sqrt(1-(cor3x+30).^2)+50;

cor4x = 29:0.1:30;
cor4y = -sqrt(1-(cor4x-30).^2)+50;

outside_x = [30 30 -30 -30 30];
outside_y = [50 -50 -50 50 50];

goal1_x = [-3.65 -3.65 3.65 3.65];
goal1_y = [-50 -50 -50 -50];
goal1_z = [0 2.4 2.4 0];
goal2_x = [-3.65 -3.65 3.65 3.65];
goal2_y = [50 50 50 50];
goal2_z = [0 2.4 2.4 0];

goalkeeper1_x = [20.15 20.15 -20.15 -20.15];
goalkeeper1_y = [-50 -33.5 -33.5 -50];
goalkeeper2_x = [20.15 20.15 -20.15 -20.15];
goalkeeper2_y = [50 33.5 33.5 50];

goal_inner1_x = [9.15 9.15 -9.15 -9.15];
goal_inner1_y = [-50 -44.5 -44.5 -50];
goal_inner2_x = [9.15 9.15 -9.15 -9.15];
goal_inner2_y = [50 44.5 44.5 50];

midline_x = [-30 30];
midline_y = [0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Plot football field and animate the trajectory of football
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O_1 = animatedline('Color','r');
hold on;

patch([35 35 -35 -35], [55 -55 -55 55], [0 0 0 0], 'green') 

plot([30 30 -30 -30 30],[50 -50 -50 50 50], 'color', 'w', 'LineWidth',2)

plot3(0,0,0, '.','MarkerSize',18, 'color','w');
plot3(0,39,0, '.','MarkerSize',15, 'color','w');
plot3(0,-39,0, '.','MarkerSize',15, 'color','w');

plot3(goal1_x,goal1_y,goal1_z, 'color', 'k', 'LineWidth',4);
plot3(goal2_x,goal2_y,goal2_z, 'color', 'k', 'LineWidth',4);

plot(goalkeeper1_x, goalkeeper1_y, 'color', 'w', 'LineWidth',2);
plot(goalkeeper2_x, goalkeeper2_y, 'color', 'w', 'LineWidth',2);
plot(goal_inner1_x, goal_inner1_y, 'color', 'w', 'LineWidth',2);
plot(goal_inner2_x, goal_inner2_y, 'color', 'w', 'LineWidth',2);

plot(Rangex,Cy1, 'color', 'w', 'LineWidth',2);
plot(Rangex,Cy2, 'color', 'w', 'LineWidth',2);
plot(Range2x,Sy1, 'color', 'w', 'LineWidth',2);
plot(Range2x,Sy2, 'color', 'w', 'LineWidth',2);

plot(corx,cory, 'color', 'w', 'LineWidth',2);
plot(cor2x,cor2y, 'color', 'w', 'LineWidth',2);
plot(cor3x,cor3y, 'color', 'w', 'LineWidth',2);
plot(cor4x,cor4y, 'color', 'w', 'LineWidth',2);

plot(midline_x, midline_y, 'color', 'w', 'LineWidth',2)

U = plot3(x(1),y(1),z(1),'o','MarkerSize',7 ,'MarkerFaceColor','k');

grid on;
xlabel('x position') 
ylabel('y position')
zlabel('z position')

set(gca,'FontSize',20)
xlim([-35 35])
ylim([-55 55])
zlim([-0.1 30])
set(gcf,'units','normalized','outerposition',[0 0 1 1])

view([vp 20])
daspect([1 1 1]);
pause(3)




for i = 1:length(x)
    U.XData = x(i);
    U.YData = y(i);
    U.ZData = z(i);
    addpoints(O_1,x(i),y(i),z(i));
    drawnow update
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









