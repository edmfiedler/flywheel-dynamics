function y = rk4_tot(f,y,x1,x2,x3,x4,x5,x6,u,t,h)

    k1 = f(t,y,x1,x2,x3,x4,x5,x6,u);
    k2 = f(t+h/2,y+h*k1/2,x1,x2,x3,x4,x5,x6,u);
    k3 = f(t+h/2,y+h*k2/2,x1,x2,x3,x4,x5,x6,u);
    k4 = f(t+h,y+h*k3,x1,x2,x3,x4,x5,x6,u);
    y = y + (h/6)*(k1+2*k2+2*k3+k4);
end