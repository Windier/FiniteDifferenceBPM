clear all
close all

addpath('../functions')

%% Parameters

N = 350;
L = 80;
dx = L/(N-1);
dy = dx;
Ds = dx*dx;

x = -L/2:dx:L/2;
y = x;
[X,Y] = meshgrid(x,y);
[phi, rho] = cart2pol(X,Y);

n0 = 1.5078;
lambda = 640e-3;

k0 = 2*pi/lambda;
k = k0*n0;

dz = 20;
zFinal = 5000;

Nz = round(zFinal/dz);
zList = linspace(0,zFinal, Nz);

pmlWidth = L/8;
insideIndex = round(N*[pmlWidth/L, 1-pmlWidth/L]);
region = insideIndex(1):insideIndex(2);


%% Graphics
fig1 = figure(1);
fig1.Position = [336 155 500 500];
fig1.Color = 'w';
hold on
img1 = imagesc(x,y,zeros(N));
rectangle('Position',[-L/2+pmlWidth -L/2+pmlWidth L-2*pmlWidth L-2*pmlWidth],'LineWidth',2,'EdgeColor','r')
title('Intensity','Interpreter','latex')
cmap = nonlinearMap(inferno(1024),2);
colormap(cmap);
colorbar
axis image
xlim([-L/2+pmlWidth, L/2-pmlWidth])
ylim([-L/2+pmlWidth, L/2-pmlWidth])
set(gca,'FontSize',12, ...
    'TickLabelInterpreter','latex')
xlabel('$x(\mu \mathrm{m})$','interpreter', 'latex')
ylabel('$y(\mu \mathrm{m})$','interpreter', 'latex')


%% Initial Conditions
w0 = 10;
ang = [0.0, 0.0];
beamAng = 20*(pi/180);
beamDist = 0.0;

G = @(x,w0) exp(-x.^2./w0.^2);
LG = rho.*G(X,w0/2).*G(Y,w0/2).*exp(1i*phi);

F = G;

% createInitialField applied a function (first argument) as 
% F(X+sin(beamAng)*beamDist,w0).*F(Y-cos(beamAng)*beamDist,w0)
% and gives it an angle based on the 'ang' variable in degrees
U = createInitialField(F, X, Y, w0, k, beamAng, beamDist, ang);


%% Waveguide Construction
guideRadius = 0.8;
guideSeparation = 3.0;

rect = @(x,w0) double(abs(x/w0) < 1/2);
rho = @(x,y) sqrt(x.^2 + y.^2);

WG = @(a,b) K(X-a,Y-b);

waveguides = {@(z) curvedWaveguide(z,[0,0],[0,0], 0, zFinal, WG, 2.2e-3);};


%% Visualization
fig2 = figure(2);
fig2.Position =  [795 175 560 420];
fig2.Color = 'w';

% Creates a plot (isosurface or line) to visualize the waveguides
% parameters: 
% type: 1 (isosurface), 2 (lines)
% coloring: 1 (by segment), 2 (by refractive index change)
params = struct( ...
    'type', 1, ...
    'coloring', 2, ...
    'dz', dz, ...
    'n0', n0, ...
    'L', L, ...
    'N', N, ...
    'zFinal', zFinal);

visualize(waveguides, params);


%% Simulation
simParams = struct( ...
    'N', N, ...
    'L', L, ...
    'n0', n0, ...
    'lambda', lambda, ...
    'plotStep', 1);

U = FDpropagate(U, simParams, zList, waveguides,{img1});

function out = K(X,Y)
    rect = @(x,w0) double(abs(x/w0) < 1/2);
    % Elipticity shows as unbalanced x and y
    rho = @(x,y) sqrt(x.^2 + 0.4*y.^2); 

    numWG = 18; %num of waveguides
    radiusWG = 2; %radius of each waveguide
    distance = 10; % distance from the waveguide to the center (0,0)

    thetaStep = 2*pi/numWG;

    out = zeros(size(X));
    for i = 0:(numWG-1)
        x0 = distance*cos(i*thetaStep);
        y0 = distance*sin(i*thetaStep);
        out = out + rect(rho(X-x0,Y-y0), radiusWG);
    end
end






