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
dist = 7.5;
Nz = round(zFinal/dz);
zList = linspace(0,zFinal, Nz);

pmlWidth = L/8;

%% Graphics
insideIndex = round(N*[pmlWidth/L, 1-pmlWidth/L]);
region = insideIndex(1):insideIndex(2);
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
guideWaist = 0.8;
guideSeparation = 3.0;

rho = @(x,y) sqrt(x.^2 + y.^2);
K = @(X,Y) exp(-(rho(X,Y) - guideSeparation).^2./guideWaist.^2);

WG = @(a,b) K(X-a,Y-b);


curvedWaveguide = @curvedWaveguide;


L1 = 2000;
L2 = 2000;

waveguides = {@(z) curvedWaveguide(z,[0,0],[0,0], 0, L1, WG, 1.2e-3);
              @(z) curvedWaveguide(z,[0,0],[0,0], L1, L2, WG, 3.2e-3);
              @(z) curvedWaveguide(z,[0,0],[0,0], L1+L2, zFinal, WG, 0.5e-3);};


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








