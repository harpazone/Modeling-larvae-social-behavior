function [visual_angle, right_side] = calcVisualAngle(neigh_dist,neigh_ang,neigh_rel_ori,horz_size)
%description: calculates the visual angle that a neighbor occupies on the
%eye


%Input:  neigh_dist - Nx1 distance to all N neighbors
%        neigh_ang - angle in space from fish to neighbors [deg]
%        neigh_rel_ori - relative orientation of neighbros [deg]
%        horz_size - length of neighbor (must match units of distance)




%Output:  visual_angle - Nx1 visual angle of neighbors [deg]
%         left_right - Nx1 with 1's if neighbors are on the right visual
%                      field and 0 otherwise


%...........Local Variable definitions..........
% random data for testing:
% horz_size = 20;
% neigh_dist = rand(5,1)*50 + 20;
% neigh_ang = rand(5,1)*360;
% neigh_rel_ori = rand(5,1)*360;

%.................Main Function.................

% calc x,y position with respect to focal fish
neigh_x = sind(neigh_ang).*neigh_dist;
neigh_y = cosd(neigh_ang).*neigh_dist;

% find position of start and end points of the neighbros
neigh_x_body_strt = neigh_x + sind(neigh_rel_ori)*horz_size/2;
neigh_y_body_strt = neigh_y + cosd(neigh_rel_ori)*horz_size/2;

neigh_x_body_end = neigh_x - sind(neigh_rel_ori)*horz_size/2;
neigh_y_body_end = neigh_y - cosd(neigh_rel_ori)*horz_size/2;

% the positions are the vectors
% calc angle between vectors. 
visual_angle = angOfVectors([neigh_x_body_strt neigh_y_body_strt],...
    [neigh_x_body_end,neigh_y_body_end],1);

% decide if fish are to the left (0) or to the right (1)
right_side = neigh_ang>0;

% plot to confirm:
% for i = 1:size(neigh_x,1);
%     plot([neigh_x_body_strt(i), neigh_x_body_end(i)],...
%         [neigh_y_body_strt(i), neigh_y_body_end(i)]); 
%     hold on;
% end
% plot([0 0],[-horz_size/2 horz_size/2],'k');
% axis([-100 100 -100 100])
% 
% [neigh_dist,neigh_ang neigh_rel_ori visual_angle]

%............Call for local functions...........




