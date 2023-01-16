clear
clc

parameters

% Height
dp0 = [[-2.1:0.2:-0.1] [0.1:0.2:2.1]]*10^-3;     

Fdat = [];
for g = 1:length(dp0)
    dp = dp0(g);
    r_upper = [20 23 26]*10^-3;

    main

    Fdat = [Fdat Ftot]
end

%% Plot
plot(dp0,Fdat(2,:))