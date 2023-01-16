function [pos_std,cur_std] = main_PID_spin(Om)

% clear
% clc
% close all

addpath('functions')

%% Load parameters

parameters

Om = 2*pi*(Om)/60; % [rads^-1] Spin velocity
ts = 5e-5; % [s] Step size

% Helper variables
k = (1/4)*mu0*(n^2)*Aa;
ki = 4*k*i0*0.9238795/(s0^2);
ks = -4*k*(i0^2)*0.9238795/(s0^3);
L = mu0*(n^2)*Aa/(2*s0);
ku = -mu0*(n^2)*Aa*i0/(2*s0^2);

% Controller parameters
kp1 = 0.0032;
ti1 = 5e-7;
td1 = 0.1;

kp2 = 10000;
ti2 = 5;

% FPGA Gains
% 0.0032, 5e-7, 0.1 -- constrained to +-0.6
% 10000, 5 -- constrained to +-950

%% Nonlinear Equations of motion

% Mechanical System
dx =    @(t,x,vx,y,vy,z,vz,phi,ophi,psi,opsi,ia1x,ia2x,ia1y,ia2y) vx;
dvx =   @(t,x,vx,y,vy,z,vz,phi,ophi,psi,opsi,ia1x,ia2x,ia1y,ia2y) (1/m)*(5*cos(Om*t)+ki*ia1x+ki*ia2x-2*ks*x+kpr*(x-lp*sin(psi)));
dy =    @(t,y,vy,z,vz,phi,ophi,psi,opsi,x,vx,ia1x,ia2x,ia1y,ia2y) vy;
dvy =   @(t,y,vy,z,vz,phi,ophi,psi,opsi,x,vx,ia1x,ia2x,ia1y,ia2y) (1/m)*(5*sin(Om*t)+ki*ia1y+ki*ia2y-2*ks*y+kpr*(y+lp*sin(phi)*cos(psi)));
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

%% Butterworth filter

[Af,Bf,Cf,Df] = butter(2,2*pi*4000,'s');

b_sys = ss(Af,Bf,Cf,Df);
bd_sys = c2d(b_sys,ts);
Af = bd_sys.a; Bf = bd_sys.b;
Cf = bd_sys.c; Df = bd_sys.d;

%% Initialise states

% xp = 200e-6; vxp = 0;
% yp = 275e-6; vyp = 0;
% zp = 0; vzp = 0;
% phip = 0.0001; ophip = 0;
% psip = 0.0002; opsip = 0;

xp = 0; vxp = 0;
yp = 0; vyp = 0;
zp = 0; vzp = 0;
phip = 0; ophip = 0;
psip = 0; opsip = 0;

ia1xp = 0; ia1yp = 0;
ia2xp = 0; ia2yp = 0;

% State vector
Xp = [xp vxp yp vyp phip ophip psip opsip ia1xp ia1yp ia2xp ia2yp]';
Yp = [xp+l_s*sin(psip);
    yp-l_s*sin(phip)*cos(psip);
    xp-l_s*sin(psip);
    yp+l_s*sin(phip)*cos(psip);
    ia1xp;
    ia1yp;
    ia2xp;
    ia2yp];

Xf = Yp;
dXf = zeros(12,1);

u = [0 0 0 0]';
i_dp = [0 0 0 0]';
Yref = [0 0 0 0]';
Yrefp = Yref;
ix = 0; ixref = 0;
i_ia = 0; i_i_d = 0;

ep_pid = 0; ei_pid = 0;
ep_pi = 0; ei_pi = 0;

%% Simulation loop

sim_length = 0:ts:2;
Data = [];

for k = 1:length(sim_length)

    % Runge-Kutta 2nd order integration of Mechanical System
    [x, vx]     = rk4_2nd_tot(dx,dvx,xp,vxp,yp,vyp,zp,vzp,phip,ophip,psip,opsip,ia1xp,ia2xp,ia1yp,ia2yp,k*ts,ts);
    [y, vy]     = rk4_2nd_tot(dy,dvy,yp,vyp,zp,vzp,phip,ophip,psip,opsip,xp,vxp,ia1xp,ia2xp,ia1yp,ia2yp,k*ts,ts);
    [z, vz]     = rk4_2nd_tot(dz,dvz,zp,vzp,phip,ophip,psip,opsip,xp,vxp,yp,vyp,ia1xp,ia2xp,ia1yp,ia2yp,k*ts,ts);
    [phi, ophi] = rk4_2nd_tot(dphi,dophi,phip,ophip,psip,opsip,xp,vxp,yp,vyp,zp,vzp,ia1xp,ia2xp,ia1yp,ia2yp,k*ts,ts);
    [psi, opsi] = rk4_2nd_tot(dpsi,dopsi,psip,opsip,xp,vxp,yp,vyp,zp,vzp,phip,ophip,ia1xp,ia2xp,ia1yp,ia2yp,k*ts,ts);

    % Runge-Kutta integration of Electrical System
    ia1x = rk4_tot(dia1x,ia1xp,vxp,vyp,psip,opsip,phip,ophip,u,k*ts,ts);
    ia1y = rk4_tot(dia1y,ia1yp,vxp,vyp,psip,opsip,phip,ophip,u,k*ts,ts);
    ia2x = rk4_tot(dia2x,ia2xp,vxp,vyp,psip,opsip,phip,ophip,u,k*ts,ts);
    ia2y = rk4_tot(dia2y,ia2yp,vxp,vyp,psip,opsip,phip,ophip,u,k*ts,ts);

    % State vector
    X = [x vx y vy phi ophi psi opsi ia1x ia1y ia2x ia2y]';

    % Output vector
    Y = [x+l_s*sin(psi);
    y-l_s*sin(phi)*cos(psi);
    x-l_s*sin(psi);
    y+l_s*sin(phi)*cos(psi);
    ia1x;
    ia1y;
    ia2x;
    ia2y];%+[randn(4,1)*0.77e-6;randn(4,1)*0.002];

    % Butterworth filter
    [dXf,Xf,Yf] = butter_filter(Af,Bf,Cf,Df,Xf,dXf,Y);
    Y = Yf+diag([0.9e-6 0.9e-6 2.5e-6 0.9e-6 0.0025 0.0025 0.0025 0.0025])*randn(8,1);

    % Control Law
    [i_d,ep_pid,ei_pid] = PID_FPGA(kp1,ti1,td1,Yref,Y(1:4)*10^6,ep_pid,ei_pid);
%     if k >= 2000
%         i_d = i_d+[0 0 0 0]';
%     end
    % i_d = min(0.6*ones(4,1),max(-0.6*ones(4,1),i_d));

    [u,ep_pi,ei_pi] = PI_FPGA(kp2,ti2,i_d,Y(5:end),ep_pi,ei_pi);
    % u = min(950*ones(4,1),max(-950*ones(4,1),u));
    u = u*15/1000+diag([0.005 0.005 0.005 0.005])*randn(4,1);

    % Update variables
    xp = x; vxp = vx;
    yp = y; vyp = vy;
    zp = z; vzp = vz;
    phip = phi; ophip = ophi;
    psip = psi; opsip = opsi;

    ia1xp = ia1x; ia1yp = ia1y;
    ia2xp = ia2x; ia2yp = ia2y;

    i_dp = i_d;
    Yp = Y;
    Yrefp = Yref;

    % Calculate forces
    xa1 = x+l*sin(psi);
    ya1 = y-l*sin(phi)*cos(psi);
    xa2 = x-l*sin(psi);
    ya2 = y+l*sin(phi)*cos(psi);
    pos_a = [xa1 ya1 xa2 ya2]';
    
    F = ki*Y(5:end)-ks*pos_a;
    Data = [Data [Y;F;u;i_dp]];
end

pos_std = std(Data(1,2000:end))+std(Data(2,2000:end))+std(Data(3,2000:end))+std(Data(4,2000:end));
cur_std = std(Data(5,2000:end))+std(Data(6,2000:end))+std(Data(7,2000:end))+std(Data(8,2000:end));

% %% Plot
% close all
% 
% set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');
% set(0,'defaultTextInterpreter','latex');
% 
% b = 0.09; e = 0.25;
% 
% figure('Units','inches','Position',[7 6 6.693 4],'PaperPositionMode','auto');
% subplot(2,2,1)
% hold on
% plot(sim_length,Data(1,:)*10^6,'k')
% hold off
% title('$x_{s1}$')
% xlabel('Time [s]')
% ylabel('Position [$\mu$m]')
% grid on
% subplot(2,2,2)
% hold on
% plot(sim_length,Data(2,:)*10^6,'k')
% hold off
% title('$y_{s1}$')
% xlabel('Time [s]')
% ylabel('Position [$\mu$m]')
% grid on
% subplot(2,2,3)
% hold on
% plot(sim_length,Data(3,:)*10^6,'k')
% hold off
% title('$x_{s2}$')
% xlabel('Time [s]')
% ylabel('Position [$\mu$m]')
% grid on
% subplot(2,2,4)
% hold on
% plot(sim_length,Data(4,:)*10^6,'k')
% hold off
% title('$y_{s2}$')
% xlabel('Time [s]')
% ylabel('Position [$\mu$m]')
% grid on
% sgtitle('Time series plots of positions across sensor planes')
% 
% figure('Units','inches','Position',[7 6 6.693 4],'PaperPositionMode','auto');
% xc = -300e-6:300e-6/100:300e-6;
% yc = sqrt((300e-6^2)-xc.^2);
% plot(Data(1,:)*10^6,Data(2,:)*10^6)
% hold on
% plot(Data(3,:)*10^6,Data(4,:)*10^6)
% plot(xc*10^6,yc*10^6,'k--')
% plot(xc*10^6,-yc*10^6,'k--')
% hold off
% xlim([-350 350])
% ylim([-350 350])
% grid('on')
% axis('equal')
% legend('Top','Bottom','Boundary')
% xlabel('Position in x [$\mu$m]')
% ylabel('Position in y [$\mu$m]')
% title('Sensed position of top and bottom of the rotor in relation to boundary')
% 
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
% 
% figure('Units','inches','Position',[7 6 6.693 4],'PaperPositionMode','auto');
% sgtitle('AMB Control Currents in the orthogonal plane directions')
% subplot(2,2,1)
% hold on
% plot(sim_length,Data(5,:),'k')
% % plot(sim_length,Data(17,:),'k--')
% hold off
% title('$i_{a1x}$')
% xlabel('Time [s]')
% ylabel('Current [A]')
% grid on
% subplot(2,2,2)
% hold on
% plot(sim_length,Data(6,:),'k')
% % plot(sim_length,Data(18,:),'k--')
% hold off
% title('$i_{a1y}$')
% xlabel('Time [s]')
% ylabel('Current [A]')
% grid on
% subplot(2,2,3)
% hold on
% plot(sim_length,Data(7,:),'k')
% % plot(sim_length,Data(19,:),'k--')
% hold off
% title('$i_{a2x}$')
% xlabel('Time [s]')
% ylabel('Current [A]')
% grid on
% subplot(2,2,4)
% hold on
% plot(sim_length,Data(8,:),'k')
% % plot(sim_length,Data(20,:),'k--')
% hold off
% title('$i_{a2y}$')
% xlabel('Time [s]')
% ylabel('Current [A]')
% grid on
% % legend('Actual','Desired')
% 
% figure('Units','inches','Position',[7 6 6.693 4],'PaperPositionMode','auto');
% subplot(2,2,1)
% hold on
% stairs(sim_length,Data(1,:)*10^6,'k')
% hold off
% title('$x_{s1}$')
% xlabel('Time [s]')
% ylabel('Position [$\mu$m]')
% grid on
% xlim([b e])
% subplot(2,2,2)
% hold on
% stairs(sim_length,Data(3,:)*10^6,'k')
% hold off
% title('$x_{s2}$')
% xlabel('Time [s]')
% ylabel('Position [$\mu$m]')
% grid on
% xlim([b e])
% subplot(2,2,3)
% hold on
% stairs(sim_length,Data(5,:),'k')
% % plot(sim_length,Data(17,:),'k--')
% hold off
% title('$i_{a1x}$')
% xlabel('Time [s]')
% ylabel('Current [A]')
% grid on
% xlim([b e])
% subplot(2,2,4)
% hold on
% stairs(sim_length,Data(7,:),'k')
% % plot(sim_length,Data(19,:),'k--')
% hold off
% title('$i_{a2x}$')
% xlabel('Time [s]')
% ylabel('Current [A]')
% grid on
% xlim([b e])
% % legend('Experimental','Simulation','Location','NorthEast')
% sgtitle('Response in position and current in the x-direction')
% 
% figure('Units','inches','Position',[7 6 6.693 4],'PaperPositionMode','auto');
% subplot(2,2,1)
% hold on
% stairs(sim_length,Data(2,:)*10^6,'k')
% hold off
% title('$y_{s1}$')
% xlabel('Time [s]')
% ylabel('Position [$\mu$m]')
% grid on
% xlim([b e])
% subplot(2,2,2)
% hold on
% stairs(sim_length,Data(4,:)*10^6,'k')
% hold off
% title('$y_{s2}$')
% xlabel('Time [s]')
% ylabel('Position [$\mu$m]')
% grid on
% xlim([b e])
% subplot(2,2,3)
% hold on
% stairs(sim_length,Data(6,:),'k')
% % plot(sim_length,Data(18,:),'k--')
% hold off
% title('$i_{a1y}$')
% xlabel('Time [s]')
% ylabel('Current [A]')
% grid on
% xlim([b e])
% subplot(2,2,4)
% hold on
% stairs(sim_length,Data(8,:),'k')
% % plot(sim_length,Data(20,:),'k--')
% hold off
% title('$i_{a2y}$')
% xlabel('Time [s]')
% ylabel('Current [A]')
% grid on
% xlim([b e])
% % legend('Experimental','Simulation','Location','NorthEast')
% sgtitle('Response in position and current in the y-direction')

end