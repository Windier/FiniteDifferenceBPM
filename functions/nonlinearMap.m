function out = nonlinearMap(map,scalingIntensity)
%% Colormap nao linear
N = length(map);
cMap = map;
dataMax = 1.0;
dataMin = 0;
centerPoint = 0.0;
% scalingIntensity = 2;
x = 1:length(cMap);
x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
x = scalingIntensity * x/max(abs(x));
x = sign(x).* exp(abs(x));
x = x - min(x); x = x*(N-1)/max(x)+1;
newMap = interp1(x, cMap, 1:N);
out = newMap;
end

