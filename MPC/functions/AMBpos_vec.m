function [pos_a,v_a] = AMBpos_vec(x,vx,y,vy,psi,opsi,phi,ophi,l);
    % The determined forces can be used to find desired set-point currents
    % Find positions
    xa1 = x + l*sin(psi);
    ya1 = y - l*sin(phi)*cos(psi);
    xa2 = x - l*sin(psi);
    ya2 = y + l*sin(phi)*cos(psi);

    % Obtain AMB plane velocities
    dxa1 = vx + opsi*l*cos(psi);
    dya1 = vy - ophi*l*cos(phi)*cos(psi) + opsi*l*sin(phi)*sin(psi);
    dxa2 = vx - opsi*l*cos(psi);
    dya2 = vy + ophi*l*cos(phi)*cos(psi) - opsi*l*sin(phi)*sin(psi);

    pos_a = [xa1 ya1 xa2 ya2]';
    v_a = [dxa1 dya1 dxa2 dya2]';
end

