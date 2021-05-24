function [x_rand,y_rand] = dotsInAcircle(n_dots,R,offset,center)
%This function generated n random dots in a circle of radius R
%if 'offset' exists than the points are in a ring between offset and R

%Input
    %n_dots
    %R - radius of the circle
    %center - of the circle
    
    
if nargin==2
    offset = 0;
end

if nargin==3
    center = [0 0];
end

Ro = R-offset; % reduce R by the offset

t = 2*pi*rand(n_dots,1); %random angles
r = Ro*sqrt(rand(n_dots,1)); %random distances
r = r+offset;
x_rand = center(1)+r.*cos(t);
y_rand = center(2)+r.*sin(t);
end

