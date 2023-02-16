function [field,slice] = propagateField(U,z_list,img1,L,N,lambda)
dx = L/N;
Nz = length(z_list);

pad = N;
U0 = padarray(U,[pad, pad],'post');
slice = zeros(N,Nz);
[M,~] = size(U0);

A = fftshift(fft2(fftshift(U0)));
one_over_wvl = 1/lambda.^2;
[vx vy] = meshgrid(-1/(2*dx):1/(2*L):1/(2*dx)-1/(2*L));
p = sqrt(one_over_wvl - vx.^2 - vy.^2);

for i = 1:Nz
    z = z_list(i);
    H = exp(1i*2*pi*z*p);
    Uz = fftshift(ifft2(fftshift(A.*H)));
    field = Uz(1:N,1:N);
    I = abs(field);
    I = I./max(max(I));
%     slice(:,i) = I(:,N/2);
    img1.CData = I./max(max(I));
    drawnow
end
end

