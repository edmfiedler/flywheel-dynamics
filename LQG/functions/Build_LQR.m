u_max = 150; ep_max = 10; ei_max = 2.5;
pos_max = 100e-6; v_max = 0.1;
ang_max = 5e-4; om_max = 0.1;
i_max = 0.3;
fil_max = 0.7;
out_max = 100;

Q = diag([u_max u_max u_max u_max ...
    ep_max ep_max ep_max ep_max ...
    ei_max ei_max ei_max ei_max ...
    pos_max v_max pos_max v_max ...
    ang_max om_max ang_max om_max ...
    i_max i_max i_max i_max ...
    fil_max out_max fil_max out_max fil_max out_max fil_max out_max ...
    0.1 0.3 0.1 0.3 0.1 0.3 0.1 0.3].^-2);

R = diag([0.05^-2 0.05^-2 0.031^-2 0.031^-2]);

Qi = diag([1 1 1 1]*1);

Qa = zeros(44);
Qa(1:40,1:40) = Q; Qa(41:end,41:end) = Qi;

Aa = [A zeros(size(A,1),size(C(1:4,:),1));
    -ts*C(1:4,:) eye(size(C(1:4,:),1))];
Ba = [B;zeros(size(C(1:4,:),1),size(B,2))];

K = dlqr(Aa,Ba,Qa,R);
Kx = K(:,1:40);
Ki = -K(:,41:end);