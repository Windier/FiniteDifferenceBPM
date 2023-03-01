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
w0 = 3;
ang = [0.0, 0.0];
beamAng = 90*(pi/180);
beamDist = -10.0;

m = 1;
% createInitialFieldLG gives the LG field with charge m
U = createInitialFieldLG(m, X, Y, w0, k, beamAng, beamDist, ang);


%% Waveguide Construction
guideWaist = 0.8;
guideSeparation = 3.0;

rect = @(x,w0) double(abs(x/w0) < 1/2);
rho = @(x,y) sqrt(x.^2 + y.^2);

K = @(X,Y) rect(X,12*guideWaist).*rect(Y,12*guideWaist);
WG = @(a,b) K(X-a,Y-b);

waveguides = {@(z) curvedWaveguide(z,[10,0],[0,0], 0, zFinal, WG, 1.2e-3);};


%% Visualization
% fig2 = figure(2);
% fig2.Position =  [795 175 560 420];
% fig2.Color = 'w';
% 
% % Creates a plot (isosurface or line) to visualize the waveguides
% % parameters: 
% % type: 1 (isosurface), 2 (lines)
% % coloring: 1 (by segment), 2 (by refractive index change)
% params = struct( ...
%     'type', 1, ...
%     'coloring', 2, ...
%     'dz', dz, ...
%     'n0', n0, ...
%     'L', L, ...
%     'N', N, ...
%     'zFinal', zFinal);
% 
% visualize(waveguides, params);


%% Simulation
simParams = struct( ...
    'N', N, ...
    'L', L, ...
    'n0', n0, ...
    'lambda', lambda, ...
    'plotStep', 1);

U = FDpropagate(U, simParams, zList, waveguides,{img1});








