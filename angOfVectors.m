function [theta, angle, Size] = angOfVectors(vecA, vecB, normalize)
%This function gets two unit vectors (in two dimensions) and computes:
%1. the angle(in degrees) between the vectors (from 0 to pi) 
%2. the counter clockwise angle from vecA to vecB. NOTE: vecA is the
%reference vector!!! (for example should be set to [0 1] to get the angle
%from the y axis).If we want the clockwise angle - we should switch the
%vectors

%Input: vecA, vecB - can be eithr one coor each or a vector of coor (but
%       must be of the same size (vectors of the sort: nx2 where n is the
%       number of coor)
%Output: theta - the angle between the vectors (0 to pi)
%        angle - the angle from vecA to vecB (counter clockwise)
if nargin < 3
    normalize = 0; %defult is not to normalize
end

%get vector sizes
sizeA = size(vecA,1);
sizeB = size(vecB,1);

%check validity of data 
if sizeA ~= sizeB 
    if sizeA == 1
        vecA = repmat(vecA,sizeB,1);
%         disp('resize first vector to match');
    elseif sizeB == 1
        vecB = repmat(vecB,sizeA,1);
%         disp('resize second vector to match');
    else
    disp('error in data entered');
    return;
    end
end

angle = zeros(sizeA,1);
theta = zeros(sizeA,1);
Size = zeros(sizeA,1);

if ~normalize
        t_dot = sum(vecA.*vecB,2);
        %correct small deviations above 1 and below -1
        t_dot(t_dot>1) = 1;
        t_dot(t_dot<-1) = -1;
        
        %calc angle in the range 0 - 180
        theta = acos(t_dot)*180/pi;
        
        %calc angle in the range 0 - 360
        angle = mod(atan2(vecA(:,1).*vecB(:,2) - vecB(:,1).*vecA(:,2),...
            vecA(:,1).*vecB(:,1) + vecA(:,2).*vecB(:,2)),2*pi)*180/pi;
        
%     %calc angles between every coor
%     for i = 1:sizeA
%         tempA = vecA(i,1:2);
%         tempB = vecB(i,1:2);
%         
%         %calculate theta:
%         t_dot = dot(tempA, tempB);
%         t_dot = min([t_dot,1]); %if tdot is a bit over one - correct it
%         t_dot = max([t_dot,-1]); %or if tdot is abit under -1
%         theta(i) = acos(t_dot)*180/pi;
%         
%         %calculate angle
%         angle(i) = mod(atan2(tempA(1)*tempB(2) - tempB(1)*tempA(2),...
%             tempA(1)*tempB(1) + tempA(2)*tempB(2)),2*pi)*180/pi;
%         %     if angle(i)==0
%         %         pause()
%         %     end
%         
%     end
    
else
%     disp('normalizing all vecs to length 1');
        %first normlize vectors to have length 1
        sizeA = (calculateNorm(vecA)+eps); %calculate sizes
        sizeB = (calculateNorm(vecB)+eps); 
        vecA = vecA./[sizeA sizeA]; %normlize
        vecB = vecB./[sizeB sizeB];

        %calc angles between every coor
        t_dot = sum(vecA.*vecB,2);
        %correct small deviations above 1 and below -1
        t_dot(t_dot>1) = 1;
        t_dot(t_dot<-1) = -1;
        
        %calc angle in the range 0 - 180
        theta = acos(t_dot)*180/pi;
        
        %calc angle in the range 0 - 360
        angle = mod(atan2(vecA(:,1).*vecB(:,2) - vecB(:,1).*vecA(:,2),...
            vecA(:,1).*vecB(:,1) + vecA(:,2).*vecB(:,2)),2*pi)*180/pi;
    
        Size = [sizeA sizeB];
%     for i = 1:sizeA
%         tempA = vecA(i,1:2);
%         tempB = vecB(i,1:2);
%         As = (norm(tempA)+eps);
%         Bs = (norm(tempB)+eps);
%         tempA = tempA/As;
%         tempB = tempB/Bs;
%         Size(i,1) = As; Size(i,2) = Bs;
%         
%         %calculate theta:
%         t_dot = dot(tempA, tempB);
%         t_dot = min([t_dot,1]); %if tdot is a bit over one - correct it
%         t_dot = max([t_dot,-1]); %or if tdot is abit under -1
%         theta(i) = acos(t_dot)*180/pi;
%         
%         %calculate angle
%         angle(i) = mod(atan2(tempA(1)*tempB(2) - tempB(1)*tempA(2),...
%             tempA(1)*tempB(1) + tempA(2)*tempB(2)),2*pi)*180/pi;
%         %     if angle(i)==0
%         %         pause()
%         %     end
%         
%     end
    
end