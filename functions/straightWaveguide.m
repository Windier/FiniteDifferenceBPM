 function out = straightWaveguide(z,ri,rf,zi,h,fun,deltaN)

    rdif = rf - ri;
    
    if and((z>=zi),(z<zi+h))
        center = ri + rdif.*(z-zi)/(zi+h);
        out = struct('active',true, 'center', center, 'fun', fun, 'zlim', [zi, zi+h], 'deltaN', deltaN);
        return
    else 
        center = [NaN, NaN];
        out = struct('active', false, 'center', center, 'fun', fun, 'zlim', [zi, zi+h], 'deltaN', deltaN);
        return
    end
end


