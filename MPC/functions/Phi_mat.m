function Phi = Phi_mat(A,C,n)  
    Phi = zeros(n*size(C,1),size(A,2));
    Phi_i = C*A;
    for i = 1:size(C,1):n*size(C,1)
        Phi(i:i+size(C,1)-1,:) = Phi_i;
        Phi_i = Phi_i*A;
    end
end

