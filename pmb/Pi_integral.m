function Pi = Pi_integral(n_2,k_2)
    %PI_INTEGRAL Evaulation of Pi(n_2,k_2) integral
    
    %fun = @(th,n,k) 1./((1-n.*sin(th).^2).*sqrt(1-k.*sin(th).^2));
    %Pi = integral(@(th) fun(th,n_2,k_2),0,pi/2);

    Pi = ellipticPi(n_2,k_2);
end

