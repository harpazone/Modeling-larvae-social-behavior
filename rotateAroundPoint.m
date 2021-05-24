function rotated_data = rotateAroundPoint(temp_points, center, angle)
%This function rotates 2D data counter clockwise around some given point
%as a center (data can be rotated around the origin if the center is set to
%be [0 0]
%Input: temp_points - a matrix of Nx2 coor
%       center - a center point to rotate around
%       angle - an angle in radians to rotate the set counterclock wise

%Output: rotated_data - an Nx2 matrix of the rotated data

%make sure matrix is 
if size(temp_points,2) > size(temp_points,1)
    temp_points = temp_points';
end

% create a matrix from the center point which will be used later in calculations
center = repmat(center,size(temp_points,1),1);

% define a counter-clockwise rotation matrix with the given angle
R = [cos(angle) -sin(angle); sin(angle) cos(angle)];

% do the rotation...
s = temp_points - center; % shift points in the plane so that the center of rotation is at the origin
so = R*s'; % apply the rotation about the origin
rotated_data = so' + center; % shift again so the origin goes back to the desired center of rotation

% make a plot
% plot(temp_points(:,1),temp_points(:,2),'k*', rotated_data(:,1), rotated_data(:,2), 'r.',...
%     center(1,1), center(1,2), 'bo');
% axis equal