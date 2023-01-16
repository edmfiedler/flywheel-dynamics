function Bz = Bz(r,z,M,gam1,a1,K1,rp,Pi1,gam0,a0,K0,Pi0)
    % Axial magnetic flux Bz(r,z)

    Bz = -2*(10^-7)*M*((gam1/a1)*(K1+((rp-r)/(rp+r))*Pi1)-(gam0/a0)*(K0+((rp-r)/(rp+r))*Pi0));
end

