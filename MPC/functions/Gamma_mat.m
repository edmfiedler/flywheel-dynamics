function Gamma = Gamma_mat(A,B,C,n)
    
    Gamma = zeros(n*size(C,1),n*size(B,2));
    Gamma_vec = markovvec(A,B,C,n);

    for i = 0:1:(n-1)
        zs = zeros(i*size(C,1),size(B,2));
        Gamma(:,i*size(B,2)+1:i*size(B,2)+size(B,2)) = [zs;Gamma_vec(1:end-i*size(C,1),:)];
    end
end

