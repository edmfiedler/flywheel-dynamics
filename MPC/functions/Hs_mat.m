function Hs = Hs_mat(S,n)
    A = eye(n);
    B = eye(n-1);
    diag_2s = kron(A,2*S);
    diag_1s = kron(B,-S);

    z_mat = zeros(size(diag_2s));
    pad1 = z_mat;
    pad1(1:end-size(S,1),size(S,2)+1:end) = diag_1s;

    pad2 = z_mat;
    pad2(size(S,1)+1:end,1:end-size(S,1)) = diag_1s;

    sub_s = z_mat;
    sub_s(end-size(S,1)+1:end,end-size(S,2)+1:end) = S;

    Hs = diag_2s+pad1+pad2-sub_s;
end

