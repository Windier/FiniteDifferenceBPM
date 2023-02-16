function out = visualize(waveguides, parameters)

% Function to visualize the waveguides passed as a single cell array

dz = parameters.dz;
type = parameters.type;
n0 = parameters.n0;

N = parameters.N;
L = parameters.L;
zFinal = parameters.zFinal;
d = L/8;

dx = L/(N-1);
x = -L/2:dx:L/2;
y = x;

numWaveguides = length(waveguides);
colors = linspecer(numWaveguides);
for id = 1:numWaveguides

    waveguide = waveguides{id};
    zmin = waveguide(0).zlim(1);
    zmax = waveguide(0).zlim(2);
    Nz = floor((zmax - zmin)/dz);

    z = linspace(zmin, zmax, Nz);

    switch type
        case 1
            if ~isfield(parameters, 'level')
                error('Type 1 (isosurface) requires the field "level" for the isosurface level');
            end
            nVol = zeros(N, Nz, N, 'single');
            [xVol, zVol, yVol] = meshgrid(z,x,y);
            for i = 1:Nz
                % waveguides is passed as waveguides(id) instead of
                % waveguides{id} because (id) keeps it as a cell array.
                % #MATLABthings
                nVol(:,i,:) = applyRefIndex(zmin + (i-1)*dz, waveguides(id), n0);
            end
            [faces, verts] = isosurface(xVol, yVol, zVol, nVol, parameters.level);
            p = patch('Vertices',verts,'Faces',faces);
        
            % Set the face color
            set(p, 'FaceColor', colors(id,:));
        
            % Set other properties as desired, e.g. edge color, transparency
            set(p, 'EdgeColor', 'none');
            set(p, 'FaceAlpha', 0.7);
            set(p, 'FaceLighting', 'gouraud');
        case 2
            x = zeros(1,Nz);
            y = zeros(1,Nz);
            for i = 1:Nz
                center = waveguides{id}(zmin + (i-1)*dz).center;
                x(i) = center(1);
                y(i) = center(2);
            end

            hold on
            plot3(z,x,y,'-', 'LineWidth', 3, 'Color',colors(id,:))
            hold off

        otherwise
            error('Unrecognized type')

    end

end

% Add lighting and set material properties
camlight
material([0.8 0.3 0.9])

% Set view and camera projection
view(40,20)
camproj('perspective')

% Set other graph properties
xyTicks = round(linspace(-L/2+d, L/2-d, 5));
xyLim = round([-L/2+d, L/2-d]);
ax = gca;
set(ax,'XGrid','on', 'YGrid','on', 'ZGrid','on','Box', 'on', 'BoxStyle', 'full',...
    'XLim', [0,zFinal], 'YLim', xyLim,'ZLim', xyLim, ...
    'TickLabelInterpreter', 'latex', 'ZTick',xyTicks, 'YTick',xyTicks, 'FontSize', 12)
axis vis3d
axis square
pbaspect([3,1,1])
xlabel('$z\ (\mathrm{\mu m})$','Interpreter','latex')
ylabel('$x\ (\mathrm{\mu m})$','Interpreter','latex')
zlabel('$y\ (\mathrm{\mu m})$','Interpreter','latex')

end