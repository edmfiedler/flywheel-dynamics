clear
clc

parameters

% Height
z0 = [0.01:0.5:15.01]*10^-3;     

Fdat = [];
for gg = 1:length(z0)
    z_c = z0(gg);
    z_h = [z_c z_c+h/2 z_c+h];

    main

    Fdat = [Fdat Ftot]
end
