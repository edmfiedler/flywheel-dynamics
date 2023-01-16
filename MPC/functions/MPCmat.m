function [H,Mx0,Mr,Mu_1,Md,Aqp,Phi,Gamma_d] = MPCmat(A,B,C,Bv,Q,S,n)
    
    % Obtain Phi Matrix
    Phi = Phi_mat(A,C,n);

    % Obtain Gamma Matrix
    Gamma = Gamma_mat(A,B,C,n);

    % Obtain Gamma_d Matrix
    Gamma_d = Gamma_mat(A,Bv,C,n);

    % Obtain Qz matrix
    Qz = Qz_mat(Q,n);

    % Obtain Hs matrix
    Hs = Hs_mat(S,n);

    % Obtain Lam matrix
    Lam = Lam_mat(size(B,2),n);

    % Compute matrices
    H = Gamma'*Qz*Gamma+Hs;
    Mx0 = Gamma'*Qz*Phi;
    Mr = -Gamma'*Qz;
    Md = Gamma'*Qz*Gamma_d;
    Mu_1 = -[S;zeros(size(S,1)*n-size(S,1),size(S,2))];
    Aqp = [Lam;Gamma];

end

