%% Auxiliary Functions
function [U] = createInitialFieldLG(m, X, Y, w0, k, beamAng, beamDist, ang)

G = @(x,w0) exp(-x.^2./w0.^2);
Ynew = Y-cos(beamAng)*beamDist;
Xnew = X+sin(beamAng)*beamDist;
[phi, rho] = cart2pol(Xnew,Ynew);

U = (rho./w0).^(abs(m)).*G(X+sin(beamAng)*beamDist,w0).*...
    G(Y-cos(beamAng)*beamDist,w0).*...
    exp(-(tan(ang(1)*pi/180)*k*1i*X)).*...
    exp(-(tan(ang(2)*pi/180)*k*1i*Y)).*exp(1i*m*phi);

% Ensure zero boundaries
U(:,1) = 0;
U(:,end) = 0;
U(1,:) = 0;
U(end,:) = 0;

end