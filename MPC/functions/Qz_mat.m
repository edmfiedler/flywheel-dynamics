function Qz = Qz_mat(Q,n)
    A = eye(n);
    Qz = kron(A,Q);
end

