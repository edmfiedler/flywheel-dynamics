function [u,e,ei] = PI_FPGA(P,I,sp,pv,ep,ei)
    e = sp-pv;
    ei = ei + (e+ep)/2;

    u = P*e+I*ei;
end

