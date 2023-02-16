function [modeSystem,eigenValues] = findModes(numModes, n, n0, dx, k0, deltaN)

    Ds = dx*dx;
    
    [~,N] = size(n);
    M = N*N;
    id1 = (1:(M-1))';
    id2 = (1:(M+1))';

    maskld = (1 - heaviside(mod(id1-1,N)-(N-1.5)));
    maskud = (1 - heaviside(mod(id2-2,N)-(N-1.5)));
    
    B = k0.^2*(n.^2 - (n0+deltaN).^2);
    bflat = reshape(B,[M,1]);
    
    dz = 1e-4;
    r = dz/(2i*k0*(n0+deltaN));

    ax = 1.0001*r/(dx*dx);
    ay = 1.0000*r/(dx*dx);
    
    md = spdiags(-repmat([ax 2*ax*ones(1,N-2) ax],1,N).' ...
            - [ay*ones(1,N), 2*ay*ones(1,N*(N-2)), ay*ones(1,N)].' ...
            + r*bflat, 0, M, M);
%     md = spdiags(-r*(4/Ds)*ones(M,1) + r*bflat, 0, M, M);
    ldx = spdiags(ax*maskld.*ones(M-1,1), -1, M, M);
    udx = spdiags(ax*maskud.*ones(M+1,1), 1, M, M);
    ldy = spdiags(ay*ones(M-N,1), -N, M, M);
    udy = spdiags(ay*ones(M,1), N, M, M);
    
    [modeSystem,eigenValues] = eigs(md + ldx + udx + ldy + udy, numModes, 1);

end