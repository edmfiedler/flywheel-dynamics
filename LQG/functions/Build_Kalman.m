Qk = diag(1*[10 10 10 10]);
Rk = diag(1*[0.9 0.9 2.5 0.9 0.1 0.1 0.1 0.1]);
[~,Kf,~] = idare(A',C',B*Qk*B',Rk);
Kf = Kf';