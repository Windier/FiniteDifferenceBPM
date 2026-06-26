function [U] = FDpropagate(U0, parameters, zList, waveguide, varargin)

U0 = complex(single(U0));

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

skip = false;
if isfield(parameters, 'skipProp')
    skip = parameters.skipProp;
end

makeVideo = false;
if isfield(parameters, 'makeVideo')
    makeVideo = parameters.makeVideo;
    if ~isfield(parameters, 'filename')
        warning('No video filename given');
        makeVideo = false;
    end
    if ~isfield(parameters, 'fig')
        warning('No figure provided');
        makeVideo = false;
    end
end

rect = @(x,w0) double(abs(x/w0) < 1/2);
% Hminus = @(x) 0.5*(1 - sign(x)).*x;
% Hplus = @(x) 0.5*(1 + sign(x)).*x;

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

pmlOrder = 2;
sigmax = (rect(X-L/2,2*d).*sigma0x.*((X-L/2+d)/d).^pmlOrder + ...
    rect(X+L/2,2*d).*sigma0x.*((-X-L/2+d)/d).^pmlOrder);
sigmay = (rect(Y-L/2,2*d).*sigma0y.*((Y-L/2+d)/d).^pmlOrder + ...
    rect(Y+L/2,2*d).*sigma0y.*((-Y-L/2+d)/d).^pmlOrder);

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

Txplus = single(R.*Txplus(:));
Txminus = single(R.*Txminus(:));
Rx = single(R.*Rx(:));
Typlus = single(R.*Typlus(:));
Tyminus = single(R.*Tyminus(:));
Ry = single(R.*Ry(:));

Bi1 = zeros(M,1,'single');
BiN = zeros(M,1,'single');
B1j = zeros(M,1,'single');
BNj = zeros(M,1,'single');
I0 = ones(M2,1,'single');
I = ones(M2,1,'single');

ldFirst = single(-Txminus.*maskld);
udFirst = single(-Txplus.*maskud);
ldSecond = single(-Tyminus2.*maskld);
udSecond = single(-Typlus2.*maskud);

tyminus = vertcat(Bi1,Tyminus(1:end-M));
typlus =  vertcat(Typlus(M+1:end),BiN);
txminus = vertcat(Bi1,Txminus2(1:end-M));
txplus =  vertcat(Txplus2(M+1:end),BiN);

id = 1:M2;
idminus = mod(id-M-1,M2)+1;
idplus = mod(id+M-1,M2)+1;


U = U0;

t1Start = tic;

Nz = length(zList);

time = 10;

if makeVideo
    filename = sprintf('%s.mp4',parameters.filename);
    video = VideoWriter(filename,'MPEG-4');
    video.Quality = 100;
    video.FrameRate = round((Nz/plotStep)/time);
    open(video)
end
for m = 1:Nz
    z = zList(m);
    time_n0 = tic;
    n = applyRefIndex(z, waveguide, n0);
    time_n0_end = toc(time_n0);

    if mod(m, plotStep) == 0
        if plotResults

            switch numImages
                case 1
                    imageHandlers{1}.CData = abs(U).^2;
                case 2
                    imageHandlers{1}.CData = abs(U).^2;
                    imageHandlers{2}.CData = angle(U);
                case 3
                    imageHandlers{1}.CData = abs(U).^2;
                    imageHandlers{2}.CData = angle(U);
                    imageHandlers{2}.CData = n;

                otherwise

            end

            drawnow
        end
        if makeVideo
            F = getframe(parameters.fig);
            writeVideo(video,F);
        end
    end
    time_sim = tic;
    if skip
        continue
    end
    % Find the approximate values at the half timestep
    %% First Half Step

    % Boundary Condition
    % at step m
    % Bi1(:) = 0*U(2:N-1,1);
    % BiN(:) = 0*U(2:N-1,N);



    B = 0.5*k0.^2*(n.^2 - n0.^2);
    bflat = reshape(B(2:N-1,2:N-1),[M2,1]);

    GammaPlus  = complex((I + Rx - Rn*bflat));
    GammaMinus = complex((I - Ry + Rn*bflat));
    md = GammaPlus;

    u = U(2:N-1,2:N-1);
    v = complex(u(:));

%     b = tyminus.*v(idminus) + GammaMinus.*v + typlus.*v(idplus);
%     U(2:N-1,2:N-1) = reshape(thomas_fast(ldFirst,md,udFirst, ...
%         b, ...
%         M),[M,M]);

    res = thomas( ...
        ldFirst, ...
        md, ...
        udFirst, ...
        tyminus, ...
        GammaMinus, ...
        typlus, ...
        v, ...
        M);

    U(2:N-1,2:N-1) = reshape( res ,[M,M]);


    %% Second Half Step

    % Boundary Condition
    % at step m + 1/2
    % B1j(:) = U(1,2:N-1);
    % BNj(:) = U(N,2:N-1);



    bflat = reshape(transpose(B(2:N-1,2:N-1)),[M2,1]);
    GammaPlus  =   complex((I + Ry2 - Rn*bflat));
    GammaMinus =   complex((I - Rx2 + Rn*bflat));
    md = GammaPlus;

    u = transpose(U(2:N-1,2:N-1));
    v = complex(u(:));

%     b = txminus.*v(idminus) + GammaMinus.*v + txplus.*v(idplus);
%     U(2:N-1,2:N-1) = transpose(reshape(thomas_fast(ldSecond,md,udSecond, b, M),[M,M]));

    res = thomas( ...
        ldFirst, ...
        md, ...
        udFirst, ...
        txminus, ...
        GammaMinus, ...
        txplus, ...
        v, ...
        M);
    
    U(2:N-1,2:N-1) = transpose(reshape( res ,[M,M]));

    time_sim_end = toc(time_sim);
    fprintf('z: %.1f micrometer, n calc: %d ms, sim: %d ms\n', m*dz, round(1000*time_n0_end), round(1000*time_sim_end));
    % if mod(floor(100*100*m/Nz)/100, 10) == 0
        % fprintf('%.2f %%\n', 100*m/Nz);
        % fprintf('z: %.1f micrometer, n calc: %d ms, sim: %d ms\n', m*dz, round(1000*time_n0_end), round(1000*time_sim_end));
    % end

end
if makeVideo
    close(video)
end

t1End = toc(t1Start);

end