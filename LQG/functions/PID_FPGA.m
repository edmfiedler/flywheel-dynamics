function [u,e,ei] = PID_FPGA(P,I,D,sp,pv,ep,ei)
    e = sp-pv;
    ei = ei + (e+ep)/2;
    ed = e - ep;
    
    u = P*e+I*ei+D*ed;
end

