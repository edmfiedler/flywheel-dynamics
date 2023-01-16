function [x1,x2] = rk4_2nd_tot(f1,f2,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,t,h)

    k1x1 = f1(t,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14);
    k1x2 = f2(t,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14);
    k2x1 = f1(t+h/2,x1+h*k1x1/2,x2+h*k1x2/2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14);
    k2x2 = f2(t+h/2,x1+h*k1x1/2,x2+h*k1x2/2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14);
    k3x1 = f1(t+h/2,x1+h*k2x1/2,x2+h*k2x2/2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14);
    k3x2 = f2(t+h/2,x1+h*k2x1/2,x2+h*k2x2/2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14);
    k4x1 = f1(t+h,x1+h*k3x1,x2+h*k3x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14);
    k4x2 = f2(t+h,x1+h*k3x1,x2+h*k3x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14);
    
    x1 = x1 + (h/6)*(k1x1+2*k2x1+2*k3x1+k4x1);
    x2 = x2 + (h/6)*(k1x2+2*k2x2+2*k3x2+k4x2);
end

