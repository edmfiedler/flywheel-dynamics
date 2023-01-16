clear
clc

load('dat.mat')

figure('Units','inches','Position',[0 0 6.693 4/1],'PaperPositionMode','auto');
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

hold on
quiver(dat(:,1)*1e3,dat(:,2)*1e3,dat(:,3)*1e3,dat(:,4)*1e3,'r')
plot([0.02 0.026]*1e3,[1.5e-3 1.5e-3]*1e3,'k')
plot([0.02 0.026]*1e3,[-1.5e-3 -1.5e-3]*1e3,'k')
plot([0.02 0.02]*1e3,[1.5e-3 -1.5e-3]*1e3,'k')
plot([0.026 0.026]*1e3,[1.5e-3 -1.5e-3]*1e3,'k')
plot([0.023 0.023]*1e3,[1.5e-3 -1.5e-3]*1e3,'k')
hold off
axis equal
ylim([-3e-3 3e-3]*1e3)
xlim([0.018 0.027]*1e3)

xlabel('Radius from centre [mm]')
ylabel('Height from centre of magnet [mm]')

sgtitle('2-dimensional magnetic field around PMB stator')