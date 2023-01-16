clear
clc
close all

load('Fdat.mat')

% Height
dp0 = [[-2.1:0.2:-0.1] [0.1:0.2:2.1]]*10^-3;  

P = polyfit(dp0,Fdat(2,:),10);

x_p = [-2:0.1:2]*10^-3;
y_p = polyval(P,x_p);

%% Linearised function

s0 = 0;

deriv = [];
for i = 1:length(P)-1
    k = length(P)-i;
    deriv = [deriv k*P(i)*s0^(k-1)];
end

m = sum(deriv);
c = polyval(P,0)-sum(deriv)*0;

x_l = [-2:0.1:2]*1e-3;
y_l = m.*x_l+c;

%%

figure('Units','inches','Position',[0 0 6.693 4/1.2],'PaperPositionMode','auto');
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

plot(dp0*1e3,Fdat(2,:),'bx')
hold on
plot(x_p*1e3,y_p,'k')
plot(x_l*1e3,y_l,'r')
grid on
hold off

xlabel('Displacement in rotor misalignment [mm]')
ylabel('Radial force on rotor [N]')

xlim([-2e-3 2e-3]*1e3)

legend('Calculated values','Interpolated function','Linearised around 0 m','Location','SouthEast')

sgtitle('Change in radial force with varying misalignment')