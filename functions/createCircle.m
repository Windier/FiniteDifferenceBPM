function [x,y] = createCircle(r,c)
theta = linspace(0,2*pi,256);
x = c(1) + r*cos(theta);
y = c(2) + r*sin(theta);
end


