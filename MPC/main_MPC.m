clear
clc
close all

addpath('functions')

%% Load parameters

parameters

Om = 2*pi*(2000)/60; % [rads^-1] Spin velocity
ts = 5e-5; % [s] Step size

% Helper variables
k = (1/4)*mu0*(n^2)*Aa;
ki = 4*k*i0*0.9238795/(s0^2);
ks = -4*k*(i0^2)*0.9238795/(s0^3);
L = mu0*(n^2)*Aa/(2*s0);
ku = -mu0*(n^2)*Aa*i0/(2*s0^2);

%% Nonlinear Equations of motion

% Mechanical System
dx =    @(t,x,vx,y,vy,z,vz,phi,ophi,psi,opsi,ia1x,ia2x,ia1y,ia2y) vx;
dvx =   @(t,x,vx,y,vy,z,vz,phi,ophi,psi,opsi,ia1x,ia2x,ia1y,ia2y) (1/m)*(ki*ia1x+ki*ia2x-2*ks*x+kpr*(x-lp*sin(psi)));
dy =    @(t,y,vy,z,vz,phi,ophi,psi,opsi,x,vx,ia1x,ia2x,ia1y,ia2y) vy;
dvy =   @(t,y,vy,z,vz,phi,ophi,psi,opsi,x,vx,ia1x,ia2x,ia1y,ia2y) (1/m)*(ki*ia1y+ki*ia2y-2*ks*y+kpr*(y+lp*sin(phi)*cos(psi)));
dz =    @(t,z,vz,phi,ophi,psi,opsi,x,vx,y,vy,ia1x,ia2x,ia1y,ia2y) vz;
dvz =   @(t,z,vz,phi,ophi,psi,opsi,x,vx,y,vy,ia1x,ia2x,ia1y,ia2y) (1/m)*(kpa*(z-lp*cos(phi)*cos(psi)+lp));
dphi =  @(t,phi,ophi,psi,opsi,x,vx,y,vy,z,vz,ia1x,ia2x,ia1y,ia2y) ophi;
dophi = @(t,phi,ophi,psi,opsi,x,vx,y,vy,z,vz,ia1x,ia2x,ia1y,ia2y) (1/I)*(-Izz*Om*opsi+l*(-ki*ia1y+ki*ia2y-2*ks*l*sin(phi)*cos(psi))+kpr*lp*(y+lp*sin(phi)*cos(psi)));
dpsi =  @(t,psi,opsi,x,vx,y,vy,z,vz,phi,ophi,ia1x,ia2x,ia1y,ia2y) opsi;
dopsi = @(t,psi,opsi,x,vx,y,vy,z,vz,phi,ophi,ia1x,ia2x,ia1y,ia2y) (1/I)*(Izz*Om*ophi+l*(ki*ia1x-ki*ia2x-2*ks*l*sin(psi))-kpr*lp*(x-lp*sin(psi)));

% Electrical System
dia1x = @(t,ia1x,vx,vy,psi,opsi,phi,ophi,u) (1/L)*(u(1)-ia1x*R-ku*(vx+opsi*l*cos(psi)));
dia1y = @(t,ia1y,vx,vy,psi,opsi,phi,ophi,u) (1/L)*(u(2)-ia1y*R-ku*(vy-ophi*l*cos(phi)*cos(psi)+opsi*l*sin(phi)*sin(psi)));
dia2x = @(t,ia2x,vx,vy,psi,opsi,phi,ophi,u) (1/L)*(u(3)-ia2x*R-ku*(vx-opsi*l*cos(psi)));
dia2y = @(t,ia2y,vx,vy,psi,opsi,phi,ophi,u) (1/L)*(u(4)-ia2y*R-ku*(vy+ophi*l*cos(phi)*cos(psi)-opsi*l*sin(phi)*sin(psi)));

%% Build model for control and estimation

FESS_state_space

A = sysd.a; B = sysd.b; C = sysd.c; D = sysd.d;

%% Kalman filter

Qk = diag([0.1 0.1 0.1 0.1]);
Rk = diag([0.9 0.9 2.5 0.9 0.1 0.1 0.1 0.1]);
[~,Kf,~] = idare(A',C',B*Qk*B',Rk);
Kf = Kf';

%% Compute MPC related matrices

% Depth
n = 100;

% Weighting matrices
Q = [1 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0;
    0 0 0 1 0 0 0 0;
    0 0 0 0 0.01 0 0 0;
    0 0 0 0 0 0.01 0 0;
    0 0 0 0 0 0 0.01 0;
    0 0 0 0 0 0 0 0.01]*1e11;

R = [1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];

% Build Matrices for Control Law
[H,Mx0,Mr,Mu_1,Md,Aqp,Phi,Gamma_d] = MPCmat(A,B,C,B,Q,R,n);
H = (H+H')/2;

ref = 0*ones(n*size(C,1),1);

% Input upper and lower bounds
U_min = repmat(-15*ones(4,1),n,1);
U_max = -U_min;

%% Initialise states

xp = 1e-5; vxp = 0.003;
yp = -1e-5; vyp = -0.003;
phip = 1e-4; ophip = 0.022;
psip = 1e-4; opsip = 0.022;

ia1xp = 0.006; ia1yp = 0.006;
ia2xp = 0.0024; ia2yp = 0.0024;

xp = [xp vxp yp vyp phip ophip psip opsip ia1xp ia1yp ia2xp ia2yp]';
xh = xp;

% State vector
% Yp = [xp+l_s*sin(psip);
%     yp-l_s*sin(phip)*cos(psip);
%     xp-l_s*sin(psip);
%     yp+l_s*sin(phip)*cos(psip);
%     ia1xp;
%     ia1yp;
%     ia2xp;
%     ia2yp];

u = [0 0 0 0]';
i_dp = [0 0 0 0]';

%% Simulation loop

sim_length = 0:ts:0.02;
Data = [];

for k = 1:length(sim_length)

%     % Runge-Kutta 2nd order integration of Mechanical System
%     [x, vx]     = rk4_2nd_tot(dx,dvx,xp,vxp,yp,vyp,zp,vzp,phip,ophip,psip,opsip,ia1xp,ia2xp,ia1yp,ia2yp,k*ts,ts);
%     [y, vy]     = rk4_2nd_tot(dy,dvy,yp,vyp,zp,vzp,phip,ophip,psip,opsip,xp,vxp,ia1xp,ia2xp,ia1yp,ia2yp,k*ts,ts);
%     [z, vz]     = rk4_2nd_tot(dz,dvz,zp,vzp,phip,ophip,psip,opsip,xp,vxp,yp,vyp,ia1xp,ia2xp,ia1yp,ia2yp,k*ts,ts);
%     [phi, ophi] = rk4_2nd_tot(dphi,dophi,phip,ophip,psip,opsip,xp,vxp,yp,vyp,zp,vzp,ia1xp,ia2xp,ia1yp,ia2yp,k*ts,ts);
%     [psi, opsi] = rk4_2nd_tot(dpsi,dopsi,psip,opsip,xp,vxp,yp,vyp,zp,vzp,phip,ophip,ia1xp,ia2xp,ia1yp,ia2yp,k*ts,ts);
% 
%     % Runge-Kutta integration of Electrical System
%     ia1x = rk4_tot(dia1x,ia1xp,vxp,vyp,psip,opsip,phip,ophip,u,k*ts,ts);
%     ia1y = rk4_tot(dia1y,ia1yp,vxp,vyp,psip,opsip,phip,ophip,u,k*ts,ts);
%     ia2x = rk4_tot(dia2x,ia2xp,vxp,vyp,psip,opsip,phip,ophip,u,k*ts,ts);
%     ia2y = rk4_tot(dia2y,ia2yp,vxp,vyp,psip,opsip,phip,ophip,u,k*ts,ts);

%     % Output vector
%     Y = [x+l_s*sin(psi);
%     y-l_s*sin(phi)*cos(psi);
%     x-l_s*sin(psi);
%     y+l_s*sin(phi)*cos(psi);
%     ia1x;
%     ia1y;
%     ia2x;
%     ia2y];

    x = A*xp+B*u;
    Y = C*xp+diag([0.9 0.9 2.5 0.9 0.0025 0.0025 0.0025 0.0025])*randn(8,1);
    xp = x;

    xh = A*xh+B*u+Kf*(Y-C*xh);

%     X = [x vx y vy phi ophi psi opsi ia1x ia1y ia2x ia2y]';
    % MPC
    gu = Mx0*xh+Mr*ref+Mu_1*u;
    [u,info] = qpsolver(H,gu,U_min,U_max,[],[],[],[]);
    u = u(1:4)+diag([0.02 0.02 0.02 0.02])*randn(4,1);
    u = min(15*ones(4,1),max(-15*ones(4,1),u));

    % Update variables
%     xp = x; vxp = vx;
%     yp = y; vyp = vy;
%     zp = z; vzp = vz;
%     phip = phi; ophip = ophi;
%     psip = psi; opsip = opsi;
% 
%     ia1xp = ia1x; ia1yp = ia1y;
%     ia2xp = ia2x; ia2yp = ia2y;
% 
%     % Calculate forces
%     xa1 = x+l*sin(psi);
%     ya1 = y-l*sin(phi)*cos(psi);
%     xa2 = x-l*sin(psi);
%     ya2 = y+l*sin(phi)*cos(psi);
%     pos_a = [xa1 ya1 xa2 ya2]';
%     
%     F = ki*Y(5:end)-ks*pos_a;
    Data = [Data [Y]];   
end

Data_OG = Data;
%% Load LQR data
load LQR_comp_Data
LQR_Data = Data(:,1:401).*[10^6;10^6;10^6;10^6;1;1;1;1];
load PID_comp_Data
PID_Data = Data(:,1:401).*[10^6;10^6;10^6;10^6;1;1;1;1];
% load MPC_Data
% ts = 5e-5;
% sim_length = 0:ts:0.02;


Data = Data_OG;
%% Plot
close all

set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

figure('Units','inches','Position',[7 6 6.693 4],'PaperPositionMode','auto');
subplot(2,2,1)
hold on
% plot(sim_length,x1_m,'r')
plot(sim_length,Data(1,:),'k')
plot(sim_length,LQR_Data(1,:),'r')
plot(sim_length,PID_Data(1,:),'b')
hold off
title('$x_{s1}$')
xlabel('Time [s]')
ylabel('Position [$\mu$m]')
grid on
subplot(2,2,2)
hold on
% plot(sim_length,y1_m,'r')
plot(sim_length,Data(2,:),'k')
plot(sim_length,LQR_Data(2,:),'r')
plot(sim_length,PID_Data(2,:),'b')
hold off
title('$y_{s1}$')
xlabel('Time [s]')
ylabel('Position [$\mu$m]')
grid on
legend('MPC','LQG','PID','Location','SouthEast')
subplot(2,2,3)
hold on
% plot(sim_length,x2_m,'r')
plot(sim_length,Data(3,:),'k')
plot(sim_length,LQR_Data(3,:),'r')
plot(sim_length,PID_Data(3,:),'b')
hold off
title('$x_{s2}$')
xlabel('Time [s]')
ylabel('Position [$\mu$m]')
grid on
subplot(2,2,4)
hold on
% plot(sim_length,y2_m,'r')
plot(sim_length,Data(4,:),'k')
plot(sim_length,LQR_Data(4,:),'r')
plot(sim_length,PID_Data(4,:),'b')
hold off
title('$y_{s2}$')
xlabel('Time [s]')
ylabel('Position [$\mu$m]')
grid on
sgtitle('Change in position across the sensor planes')

figure('Units','inches','Position',[7 6 6.693 4],'PaperPositionMode','auto');
xc = -300:300/100:300;
yc = sqrt((300^2)-xc.^2);
plot(Data(1,:),Data(2,:))
hold on
plot(Data(3,:),Data(4,:))
plot(xc,yc,'k--')
plot(xc,-yc,'k--')
hold off
xlim([-350 350])
ylim([-350 350])
grid('on')
axis('equal')
legend('Top','Bottom','Boundary')
xlabel('Position in x [$\mu$m]')
ylabel('Position in y [$\mu$m]')
title('Sensed position of top and bottom of the rotor in relation to boundary')

% figure('Units','inches','Position',[7 6 6.693 4],'PaperPositionMode','auto');
% stairs(sim_length,Data(9,:))
% hold on
% stairs(sim_length,Data(10,:))
% stairs(sim_length,Data(11,:))
% stairs(sim_length,Data(12,:))
% hold off
% legend('$F_{a1,x}$','$F_{a1,y}$','$F_{a2,x}$','$F_{a2,y}$')
% title('AMB forces applied to mechanical system')
% xlabel('Time [s]')
% ylabel('Force [N]')
% grid on
% 
% figure('Units','inches','Position',[7 6 6.693 4],'PaperPositionMode','auto');
% stairs(sim_length,Data(13,:))
% hold on
% stairs(sim_length,Data(14,:))
% stairs(sim_length,Data(15,:))
% stairs(sim_length,Data(16,:))
% hold off
% legend('$u_1$','$u_2$','$u_3$','$u_4$')
% title('Voltage input to the system')
% xlabel('Time [s]')
% ylabel('Voltage [V]')
% grid on

figure('Units','inches','Position',[7 6 6.693 4],'PaperPositionMode','auto');
subplot(2,2,1)
hold on
% plot(sim_length,ix1_m,'r')
plot(sim_length,Data(5,:),'k')
plot(sim_length,LQR_Data(5,:),'r')
plot(sim_length,PID_Data(5,:),'b')
% plot(sim_length,Data(17,:),'k--')
hold off
title('$i_{a1x}$')
xlabel('Time [s]')
ylabel('Current [A]')
grid on
subplot(2,2,2)
hold on
% plot(sim_length,iy1_m,'r')
plot(sim_length,Data(6,:),'k')
plot(sim_length,LQR_Data(6,:),'r')
plot(sim_length,PID_Data(6,:),'b')
% plot(sim_length,Data(18,:),'k--')
hold off
title('$i_{a1y}$')
xlabel('Time [s]')
ylabel('Current [A]')
grid on
subplot(2,2,3)
hold on
% plot(sim_length,ix2_m,'r')
plot(sim_length,Data(7,:),'k')
plot(sim_length,LQR_Data(7,:),'r')
plot(sim_length,PID_Data(7,:),'b')
% plot(sim_length,Data(19,:),'k--')
hold off
title('$i_{a2x}$')
xlabel('Time [s]')
ylabel('Current [A]')
grid on
subplot(2,2,4)
hold on
% plot(sim_length,iy2_m,'r')
plot(sim_length,Data(8,:),'k')
plot(sim_length,LQR_Data(8,:),'r')
plot(sim_length,PID_Data(8,:),'b')
% plot(sim_length,Data(20,:),'k--')
hold off
title('$i_{a2y}$')
xlabel('Time [s]')
ylabel('Current [A]')
grid on
sgtitle('Change in current over time')