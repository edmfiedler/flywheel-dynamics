function K = K_integral(k_2)
    %K_INTEGRAL Evaluation of the K(k^2) integral

    %fun = @(t,k) 1./sqrt((1-t.^2).*(1-k.*t.^2));
    %K = integral(@(t) fun(t,k_2),0,1);

    [K,E] = ellipke(k_2);
end