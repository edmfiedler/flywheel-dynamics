clear
clc

v = 0:100:10000;
data = [];

for i = 1:length(v)
    v(i)
    [p,c] = main_LQG_spinning(v(i));
    data = [data [p;c]];
end