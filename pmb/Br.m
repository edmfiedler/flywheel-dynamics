function Br = Br(r,M,a1,k1_p,K1,a0,k0_p,K0,E1,E0)
    % Radial magnetic flux Br(r,z)

    Br = 2*(10^-7)*M*((a1/(2*r))*((1-k1_p^2)*K1-2*E1)-(a0/(2*r))*((1-k0_p^2)*K0-2*E0));
end

