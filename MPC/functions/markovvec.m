function H = markovvec(A,B,C,n)

si = size(C*B,1);
H = zeros(si*n,size(C*B,2));

% Intermediate variable
H_i = C;
for i = 1:si:n*si
    H(i:i+si-1,:) = H_i*B;
    H_i = H_i*A;
end

end

