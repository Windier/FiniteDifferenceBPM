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

Ncmap = 256;
cmap = gray(Ncmap);

numWaveguides = length(waveguides);
colorSegments = linspecer(numWaveguides);
for id = 1:numWaveguides

    waveguide = waveguides{id};
    zlim = getfield(waveguide(0),'zlim');
    zmin = zlim(1);
    zmax = zlim(2);
    Nz = floor((zmax - zmin)/dz);

    z = linspace(zmin, zmax, Nz);

    switch type
        case 1
            nVol = zeros(N, Nz, N, 'single');
            [xVol, zVol, yVol] = meshgrid(z,x,y);
            for i = 1:Nz
                % waveguides is passed as waveguides(id) instead of
                % waveguides{id} because (id) keeps it as a cell array.
                % #MATLABthings
                nVol(:,i,:) = applyRefIndex(zmin + (i-1)*dz, waveguides(id), n0);
            end

            wg = waveguide(0);
            level = n0 + 0.5*wg.deltaN;

            switch parameters.coloring
                case 1
                    [faces, verts] = isosurface(xVol, yVol, zVol, nVol, level);
                    p = patch('Vertices', verts, 'Faces',faces,'FaceColor', 'interp');
                    % Set the face color
                    set(p, 'FaceColor', colorSegments(id,:));
                case 2
                    [faces, verts, colors] = isosurface(xVol, yVol, zVol, nVol, level, wg.deltaN*ones(size(nVol)));
                    p = patch('Vertices', verts,'Faces',faces,'FaceColor', 'interp', ...
                        'FaceVertexCData', colors);
                otherwise
                    error('Unknown coloring parameter')
            end

            % Set other properties, e.g. edge color, transparency
            set(p, 'EdgeColor', 'none');
            set(p, 'FaceAlpha', 0.9);
            set(p, 'FaceLighting', 'gouraud');

        case 2
            x = zeros(1,Nz);
            y = zeros(1,Nz);
            for i = 1:Nz
                wg = waveguides{id}(zmin + (i-1)*dz);
                center = wg.center;
                x(i) = center(1);
                y(i) = center(2);
            end

            switch parameters.coloring
                case 1
                    colors = linspecer(numWaveguides);
                case 2
                    colors = cmap(floor(linspace(1,Ncmap-20,numWaveguides)), :);
                otherwise
                    error('Unknown coloring parameter')
            end

            hold on
            plot3(z,x,y,'-', 'LineWidth', 3, 'Color',colors(id,:))
            hold off

        otherwise
            error('Unrecognized type')
    end

end

if parameters.coloring == 2
    cbar = colorbar;
    ylabel(cbar,'$\Delta n$','FontSize',16, 'interpreter','latex');
    colormap(gray)
    if parameters.type == 2
        warning('Colorbar scale does not represent actual ref. index change yet')
    end
end



% Add lighting and set material properties
camlight
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