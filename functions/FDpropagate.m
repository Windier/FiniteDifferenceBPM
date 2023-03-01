function [U] = FDpropagate(U0, parameters, zList, waveguide, varargin)

% U0: Initial field distribution,

% parameters: Struct passed to the solver which are:
% N: number of pixels in one dimension,
% L: physical length of the computational window,
% n0: base refractive index,
% lambda: operating wavelength
% Example: simParams = struct('N', N, ...
%                             'L', L, ...
%                             'n0', n0, ...
%                             'lambda', lambda);

% zList, waveguide

if nargin < 4
    error('FDpropagate requires at least 4 input arguments: U0, parameters, zList, waveguide');
end

plotResults = false;
if nargin == 5
    plotResults = true;
    imageHandlers = varargin{1};
    numImages = length(imageHandlers);
    try 
        plotStep = parameters.plotStep;
    catch
        plotStep = 1;
    end
end

rect = @(x,w0) double(abs(x/w0) < 1/2);
Hminus = @(x) 0.5*(1 - sign(x)).*x;
Hplus = @(x) 0.5*(1 + sign(x)).*x;

N = parameters.N;
L = parameters.L;

dx = L/(N-1);
dy = dx;
Ds = dx*dx;

x = -L/2:dx:L/2;
y = x;
[X,Y] = meshgrid(x,y);

n0 = parameters.n0;

dz = zList(2) - zList(1);
zFinal = zList(end);

k0 = 2*pi/parameters.lambda;
k = k0*n0;

a = 2i*k0*n0/dz;
r = 2i*k0*n0/dz;

c = 299792458*10^6; % Speed of light in micrometers
freq = k0*c; % Frequency in Hz
eps0 = 8.8541878128e-12;
sigma0x = 3*eps0*freq;
sigma0y = 3*eps0*freq;
np = n0;

% Finite-Difference-Scheme parameter
alpha = 0.5; %Scheme parameter (0.0 is fully implicit. 1.0 is fully explicit)
R = alpha/(r*Ds);
S = (1-alpha)/(r*Ds);
Rn = alpha/r;
d = L/8;

sigmax = (rect(X-L/2,2*d).*sigma0x.*((X-L/2+d)/d).^3 + ...
    rect(X+L/2,2*d).*sigma0x.*((-X-L/2+d)/d).^3);
sigmay = (rect(Y-L/2,2*d).*sigma0y.*((Y-L/2+d)/d).^3 + ...
    rect(Y+L/2,2*d).*sigma0y.*((-Y-L/2+d)/d).^3);

qy = 1./(1 - 1i*sigmax./(freq*eps0*np.^2));
qx = 1./(1 - 1i*sigmay./(freq*eps0*np.^2));

M = N-2;
M2 = M*M;
id = (1:(M*M))';
heav = @(x) 0.5*(1 + sign(x));
maskud = (1 - heav(mod(id-1,M)-(M-1.5)));
maskld = (1 - heav(mod(id-2,M)-(M-1.5)));

% Indices that exclude boundaries at i,j = 1 and N
id = 2:N-1;
jd = 2:N-1;

% PML Parameters
Txplus = 0.5*qx(id,jd).*(qx(id,jd) + qx(id+1,jd));
Txminus = 0.5*qx(id,jd).*(qx(id,jd) + qx(id-1,jd));
Rx = 0.5*qx(id,jd).*(qx(id-1,jd) + 2*qx(id,jd) + qx(id+1,jd));
Typlus = 0.5*qy(id,jd).*(qy(id,jd) + qy(id,jd+1));
Tyminus = 0.5*qy(id,jd).*(qy(id,jd) + qy(id,jd-1));
Ry = 0.5*qy(id,jd).*(qy(id,jd-1) + 2*qy(id,jd) + qy(id,jd+1));

%Flat PML parameters
Txplus2 = S.*reshape(transpose(Txplus),[M2,1]);
Txminus2 = S.*reshape(transpose(Txminus),[M2,1]);
Rx2 = S.*reshape(transpose(Rx),[M2,1]);
Typlus2 = S.*reshape(transpose(Typlus),[M2,1]);
Tyminus2 = S.*reshape(transpose(Tyminus),[M2,1]);
Ry2 = S.*reshape(transpose(Ry),[M2,1]);

Txplus = R.*Txplus(:);
Txminus = R.*Txminus(:);
Rx = R.*Rx(:);
Typlus = R.*Typlus(:);
Tyminus = R.*Tyminus(:);
Ry = R.*Ry(:);

Bi1 = zeros(M,1);
BiN = zeros(M,1);
B1j = zeros(M,1);
BNj = zeros(M,1);
I0 = ones(M2,1);
I = ones(M2,1);

ldFirst = -Txminus.*maskld;
udFirst = -Txplus.*maskud;
ldSecond = -Tyminus2.*maskld;
udSecond = -Typlus2.*maskud;

%RHS:
b = complex(zeros(M2,1));

U = U0;

t1Start = tic;

Nz = length(zList);
for m = 1:Nz
    z = zList(m);
    n = applyRefIndex(z, waveguide, n0);

    if mod(m, plotStep) == 0
        if plotResults

            switch numImages
                case 1
                    imageHandlers{1}.CData = abs(U).^2;
                case 2
                    imageHandlers{1}.CData = abs(U).^2;
                    imageHandlers{2}.CData = angle(U);
                otherwise

            end

            drawnow
        end
    end
    % Find the approximate values at the half timestep
    %% First Half Step

    % Boundary Condition
    % at step m
    % Bi1(:) = 0*U(2:N-1,1);
    % BiN(:) = 0*U(2:N-1,N);

    u = U(2:N-1,2:N-1);
    v = u(:);

    B = 0.5*k0.^2*(n.^2 - n0.^2);
    bflat = reshape(B(2:N-1,2:N-1),[M2,1]);

    GammaPlus  = complex((I + Rx - Rn*bflat));
    GammaMinus = complex((I - Ry + Rn*bflat));

    md = GammaPlus;

    % Find rhs Y of vector equation
    b = complex(vertcat(Bi1,Tyminus(1:end-M).*v(1:end-M)) + ...
        GammaMinus.*v + ...
        vertcat(Typlus(M+1:end).*v(M+1:end),BiN));
    % Solve vector equation to get points for half timestep
    U(2:N-1,2:N-1) = reshape(thomas(ldFirst,md,udFirst,b),[M,M]);
    %% Second Half Step

    % Boundary Condition
    % at step m + 1/2
    % B1j(:) = U(1,2:N-1);
    % BNj(:) = U(N,2:N-1);

    u = transpose(U(2:N-1,2:N-1));
    v = u(:);

    bflat = reshape(transpose(B(2:N-1,2:N-1)),[M2,1]);
    GammaPlus  =   complex((I + Ry2 - Rn*bflat));
    GammaMinus =   complex((I - Rx2 + Rn*bflat));

    md = GammaPlus;

    % Find rhs X of vector equation
    b = complex(vertcat(Bi1,Txminus2(1:end-M).*v(1:end-M)) + ...
        GammaMinus.*v + ...
        vertcat(Txplus2(M+1:end).*v(M+1:end),BiN));

    % Solve vector equation to get points for half timestep
    U(2:N-1,2:N-1) = transpose(reshape(thomas(ldSecond,md,udSecond,b),[M,M]));

    fprintf('z: %.1f micrometer\n', m*dz)

end

t1End = toc(t1Start);

end