 function out = curvedWaveguide(z,ri,rf,zi,h,fun,deltaN)

    rdif = 2*rf - ri;
    
    phi = atan2(rdif(2),rdif(1));
    d = sqrt(sum(rdif.^2));
    
    if and((z>=zi),(z<zi+h))
        if d == 0
            center = ri;
            out = struct('active',true, 'center', center, 'fun', fun, 'zlim', [zi, zi+h], 'deltaN', deltaN);
            return
        end
    else 
        center = [NaN, NaN];
        out = struct('active', false, 'center', center, 'fun', fun, 'zlim', [zi, zi+h], 'deltaN', deltaN);
        return
    end
   
    R = (h^2 + d^2)/(4*d);
    
    xi = ri(1);
    ax = xi+R;
    ay = zi;
    bx = d-R+xi;
    by = zi+h;
    
   	if and((z>=zi),(z<=zi+h/2))
        x = -sign(d)*sqrt(R.^2 - (z-ay).^2) + ax;
    elseif and((z>zi+h/2),(z<zi+h))
        x = sign(d)*sqrt(R.^2 - (z-by).^2) + bx;
    else
        x = NaN;
    end
    center = ([cos(phi), -sin(phi); sin(phi), cos(phi)]*([x, 0] - ri)')' + ri;

    out = struct('active', true, 'center', center, 'fun', fun, 'zlim', [zi, zi+h], 'deltaN', deltaN);
end