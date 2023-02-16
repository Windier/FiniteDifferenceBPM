clear all
close all

addpath('functions')

%% Parameters

N = 300;
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
fig1.Position = [600 100 500 500];
fig1.Color = 'w';
hold on
img1 = imagesc(x,y,zeros(N));
rectangle('Position',[-L/2+pmlWidth -L/2+pmlWidth L-2*pmlWidth L-2*pmlWidth],'LineWidth',2,'EdgeColor','r')
title('Intensity','Interpreter','latex')
cmap = nonlinearMap(inferno(1024),2);
colormap(cmap);
colorbar
% caxis([0,1])
axis image
xlim([-L/2+pmlWidth, L/2-pmlWidth])
ylim([-L/2+pmlWidth, L/2-pmlWidth])
set(gca,'FontSize',12, ...
    'TickLabelInterpreter','latex')
xlabel('$x(\mu \mathrm{m})$','interpreter', 'latex')
ylabel('$y(\mu \mathrm{m})$','interpreter', 'latex')


%% Initial Conditions
w0 = 10;
ang = [1.5, 0.0];
beamAng = 20*(pi/180);
beamDist = 16.2750;

G = @(x,w0) exp(-x.^2./w0.^2);
LG = rho.*G(X,w0/2).*G(Y,w0/2).*exp(1i*phi);

F = G;

% createInitialField applied a function (first argument) as 
% F(X+sin(beamAng)*beamDist,w0).*F(Y-cos(beamAng)*beamDist,w0)
% and gives it an angle based on the 'ang' variable in degrees
U = createInitialField(F, X, Y, w0, k, beamAng, beamDist, ang);


%% Waveguide Construction
guideWaist = 0.8;
guideSeparation1 = 3.0;
guideSeparation2 = 6.0;

rho = @(x,y) sqrt(x.^2 + y.^2);
K1 = @(X,Y) exp(-(rho(X,Y) - guideSeparation1).^2./guideWaist.^2);
K2 = @(X,Y) exp(-(rho(X,Y) - guideSeparation2).^2./guideWaist.^2);

WGType1 = @(a,b) K1(X-a,Y-b);
WGType2 = @(a,b) K2(X-a,Y-b);


L1 = 1000;
L2 = 1000;
dist = 8;
%% Waveguide Construction
waveguides = ...
    {...
    @(z) curvedWaveguide(z,[2*dist,0],[dist,0],0,L1,WGType1, 2.2e-3);...
    @(z) curvedWaveguide(z,[-2*dist,0],[-dist,0],0,L1, WGType2, 2.2e-3);...
    @(z) curvedWaveguide(z,[dist,0],[dist,0],L1,L2, WGType1, 2.2e-3);...
    @(z) curvedWaveguide(z,[-dist,0],[-dist,0],L1,L2, WGType2, 2.2e-3);...
    @(z) curvedWaveguide(z,[-dist,0],[-dist,0],L1 + L2,zFinal, WGType2, 2.2e-3);...
    @(z) straightWaveguide(z,[0,2*dist],[0,-2*dist],0,zFinal,WGType1, 1.2e-3);
    };


%% Visualization
fig2 = figure(2);
fig2.Position =  [795 175 560 420];
fig2.Color = 'w';

% Creates a plot (isosurface or line) to visualize the waveguides
params = struct( ...
    'type', 1, ...
    'level', 1.5079, ...
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
    'lambda', lambda);

U = FDpropagate(U, simParams, zList, waveguides,{img1});








