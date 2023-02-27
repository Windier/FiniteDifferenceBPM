function n = applyRefIndex(z, waveguides, n0)

numWaveguides = size(waveguides,1);
if numWaveguides ~= 0
    waveguide_at_z = waveguides{1}(z);
    n = n0*ones(size(waveguide_at_z.fun(0,0)));
    for i = 1:numWaveguides
        % Get the ith waveguide and evaluate it at z
        waveguide = waveguides{i}(z);
        if waveguide.active
            Delta_n = waveguide.deltaN;
            % Get the waveguide function
            fun = waveguide.fun;
            % Add to the refractive index
            n = n + Delta_n*fun(waveguide.center(1),waveguide.center(2));
        end
    end
else
    n = n0;
    %         n = min(n,n0+Delta_n);
end

end