function [d2neigh,a2neigh,relOri] = relativeNeighborProp(fi,xt,yt,ang)
%description: This function gets positions, and headings of fish in space
%and calculate relative positions and orientation w.r.t. a single fish in
%the group

%Input: fi - index of focal fish, xt, yt - vectors of positions in cart
%coordinates, ang - vector of angle in spcae 0 is y axis.

%Output: 



%...........Local Variable definitions..........

focalx = xt(fi);
focaly = yt(fi);
focalang = ang(fi);

% define headings
xh = sind(ang);
yh = cosd(ang);


%.................Main Function.................

% rotate world so focal fish points north
[rot_c] = rotateAroundPoint([xt yt],[focalx focaly],focalang*pi/180);
% [rot_h] = rotateAroundPoint([xh yh],[0 0],focalang*pi/180);

% quiver(xt,yt,xh,yh,1); hold on;
% quiver(focalx,focaly,xh(fi,:),yh(fi,:),1,'color',[0 0 0]); 
% quiver(rot_c(:,1),rot_c(:,2),rot_h(:,1),rot_h(:,2),1); hold on;

% calcualte distance of neighvors
vec2neighx = rot_c(:,1) - rot_c(fi,1);
vec2neighy = rot_c(:,2) - rot_c(fi,2);

% distance to neighbors
d2neigh = calculateNorm([vec2neighx vec2neighy]);

% angle to neighbor
[~,a2neigh] = angOfVectors([vec2neighx vec2neighy],[0 1],1);

% relative orintation
[~,relOri] = angOfVectors([xh yh],[xh(fi) yh(fi)],1);



% quiver(zeros(5,1),zeros(5,1),rot_h(:,1),rot_h(:,2),1); hold on;

% axis image;

%............Call for local functions...........




