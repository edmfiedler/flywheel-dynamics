%% Parameters

M = 883.31*10^3;        % Magnetisation
r_lower = [20 23 26]*10^-3;   % Lower magnet positions
nr_stator = 3;          % Lower magnet discretizations
dp = 1.86e-3;              % Displacement of upper magnet at theta = 0 in r
r_upper = [20 23 26]*10^-3;   % Upper magnet positions
nr_rotor = 3;           % Upper magnet discretizations
na = 6;                 % Discretizations around the ring
h = 3*10^-3;           % Height of the magnets
z0 = 1.6*10^-3;         %
z_h = [z0 z0+h/2 z0+h];
nz = 3;
y0 = 0;