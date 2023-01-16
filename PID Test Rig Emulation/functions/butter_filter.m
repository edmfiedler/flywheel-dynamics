function [dXf,Xf,Yf] = butter_filter(Af,Bf,Cf,Df,Xf,dXf,u);
    dXf_ = zeros(size(dXf));
    Xf_ = zeros(size(Xf));
    Yf = zeros(size(u));
    for i = 1:size(Xf,1)
        x = [dXf(i);Xf(i)];
        dx = Af*x+Bf*u(i);
        y = Cf*x+Df*u(i);

        dXf_(i) = dx(1);
        Xf_(i) = dx(2);
        Yf(i) = y;
    end

    Xf = Xf_;
    dXf = dXf_;
end

