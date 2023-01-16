function Lam = Lam_mat(s,n)
    A = eye(n);
    I = eye(s);
    Lam = kron(A,I);

    B = eye(n-1);
    nI = -eye(s);
    mI = kron(B,nI);

    zs = zeros(size(Lam));
    zs(s+1:end,1:end-s) = mI;

    Lam = Lam+zs;
end

