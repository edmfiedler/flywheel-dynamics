clear
clc
close all

load('Fdat.mat')

z0 = [0.01:0.5:15.01]*10^-3;

P = polyfit(z0,Fdat(3,:),10);

x_p = [0:0.1:14.5]*10^-3;
y_p = polyval(P,x_p);

%% Linearised function

Fw = 4.82*9.81;

polyval(P,1.6e-3)


s0 = 1.6e-3;

deriv = [];
for i = 1:length(P)-1
    k = length(P)-i
    deriv = [deriv k*P(i)*s0^(k-1)];
end

m = sum(deriv)
c = polyval(P,1.6e-3)-sum(deriv)*1.6e-3

x_l = [0:0.1:5]*1e-3;
y_l = m.*x_l+c;

%%

figure('Units','inches','Position',[0 0 6.693 4/1.2],'PaperPositionMode','auto');
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

plot(z0*1e3,Fdat(3,:),'bx')
hold on
plot(x_p*1e3,y_p,'k')
plot(x_l*1e3,y_l,'r')
grid on
hold off

xlabel('Displacement in axial separation [mm]')
ylabel('Axial force on rotor [N]')

legend('Calculated values','Interpolated function','Linearisation at 1.6 mm')

xlim([0 0.007]*1e3)
ylim([0 200])

sgtitle('Force in the z-direction at varying air gaps')

