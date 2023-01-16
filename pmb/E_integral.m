function E = E_integral(k_2)
    %E_INTEGRAL Evaluation of the E(k^2) integral
    % fun = @(t,k) sqrt(1-k_2.*t.^2)./(1-t.^2);
    % E = integral(@(t) fun(t,k_2),0,1);

    [K,E] = ellipke(k_2);
end

