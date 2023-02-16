%% Auxiliary Functions
function [U] = createInitialField(F, X, Y, w0, k, beamAng, beamDist, ang)


U = F(X+sin(beamAng)*beamDist,w0).*...
    F(Y-cos(beamAng)*beamDist,w0).*...
    exp(-(tan(ang(1)*pi/180)*k*1i*X)).*...
    exp(-(tan(ang(2)*pi/180)*k*1i*Y));

% Ensure zero boundaries
U(:,1) = 0;
U(:,end) = 0;
U(1,:) = 0;
U(end,:) = 0;

end