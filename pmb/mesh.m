clear
clc

parameters

z_mesh = [-3e-3:(6e-3/24):3e-3];
r_mesh = [17.9:0.2:28.1]*10^-3;

dat = [];

for k = 1:length(z_mesh)
    z = z_mesh(k);
    for c = 1:length(r_mesh)
        r = r_mesh(c);
        Bzk = [];
        Brk = [];
        for i = 1:3
            gam1 = z-h/2;
            a1 = sqrt(gam1^2+(r_lower(i)+r)^2);
            k1 = sqrt(4*r_lower(i)*r/(a1^2));
            k1_p = sqrt(1-k1^2);
            K1 = K_integral(k1^2);
            E1 = E_integral(k1^2);
            n1_2 = 4*r_lower(i)*r/((r_lower(i)+r)^2);
            Pi1 = Pi_integral(n1_2,k1^2);

            gam0 = z+h/2;
            a0 = sqrt(gam0^2+(r_lower(i)+r)^2);
            k0 = sqrt(4*r_lower(i)*r/(a0^2));
            k0_p = sqrt(1-k0^2);
            K0 = K_integral(k0^2);
            E0 = E_integral(k0^2);
            n2_2 = 4*r_lower(i)*r/((r_lower(i)+r)^2);
            Pi0 = Pi_integral(n2_2,k0^2);

            bz = Bz(r,z,M,gam1,a1,K1,r_lower(i),Pi1,gam0,a0,K0,Pi0);
            br = Br(r,M,a1,k1_p,K1,a0,k0_p,K0,E1,E0);

            if i == 2
                Bzk = [Bzk 2*bz];
                Brk = [Brk 2*br];
            else
                Bzk = [Bzk -bz];
                Brk = [Brk -br];
            end

        end

        Bzk = sum(Bzk);
        Brk = sum(Brk);

        dat = [dat;[r z Brk Bzk]];
    end
end

%%
hold on
quiver(dat(:,1),dat(:,2),dat(:,3),dat(:,4))
plot([0.02 0.026],[1.5e-3 1.5e-3],'k')
plot([0.02 0.026],[-1.5e-3 -1.5e-3],'k')
plot([0.02 0.02],[1.5e-3 -1.5e-3],'k')
plot([0.026 0.026],[1.5e-3 -1.5e-3],'k')
plot([0.023 0.023],[1.5e-3 -1.5e-3],'k')