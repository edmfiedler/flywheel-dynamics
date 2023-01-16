function r_abs = circ_intercept(th,r,dp)
    %CIRC_INTERCEPT for a given circle displaced at theta 0 by dp, the
    %intercept at angle theta is found.
    
    m = sin(th)/cos(th);
    
    icpt = roots([(m^2)+1 -2*dp (-(r^2)+(dp^2))]);
    
    if th <= pi
        if th < pi/2
            x_icpt = icpt(1);
        else
            x_icpt = icpt(2);
        end
    else
        if th < 3*pi/2
            x_icpt = icpt(2);
        else
            x_icpt = icpt(1);
        end
    end
    
    if th <= pi
        y_icpt = sqrt((r^2)-(x_icpt-dp)^2);
    else
        y_icpt = -sqrt((r^2)-(x_icpt-dp)^2);
    end
    
    r_abs = sqrt((x_icpt^2)+(y_icpt^2));

end