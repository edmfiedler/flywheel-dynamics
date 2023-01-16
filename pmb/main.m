% clear
% clc
% 
% parameters

%% Magnetic Field Calculations

datBk = [];

Ftot = 0;
for ang = 1:na
    th = (ang*pi/3)-(pi/3)+0.0628;
    for k = 1:nz
        z = z_h(k)+h/2;
        for c = 1:nr_rotor
            %r = circ_intercept(th,r_upper(c),dp)
            r = sqrt((r_upper(c)*cos(th))^2+(dp+r_upper(c)*sin(th))^2);
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
            
            % Calculate Force
            dl = 2*pi*r_lower(c)/6
            br_s = sum(Brk);
            bz_s = sum(Bzk);
            ang_vec = [-sin(th)*dl;cos(th)*dl;0];
            mag_vec = [cos(th)*br_s;sin(th)*br_s;bz_s];
            if c == 2
                f = -2*cross(ang_vec,mag_vec)*M*h/nz;
            else
                f = cross(ang_vec,mag_vec)*M*h/nz;
            end
            Ftot = Ftot + f;
        end
    end
end

Ftot