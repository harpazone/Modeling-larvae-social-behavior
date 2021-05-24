function [x,y,Vx,Vy,Speed, angle, State, wallD] = ...
    SimulateLarvaFishGitHub(int_type, age, loadpath, varargin)
% simjualtes Larva fish according to the interaction rules indicated.




% Input: int_type - interaction type [non_social, 7 dpf...]

%        optional input: Fs - frames per sec,
%                        N - number of fish,
%                        bout_rate - in Hz
%                        arena_rad - in cm.
%                        T - total simulation time [frames]
%                        BL - fish body length [in cm]
%                        PLOT - 1-yes(default), 0-no


% Output: x,y - position [NxT matrices]
%         Vx,Vy - velocitiy [NxT matrices]
%         speed - norm of velocity [NxT]
%         angle - heading angle [in degrees]
%         State - indicates which behavior was axecuted ....
%         WallD - distance to closes wall [in cm];


% set defaults
defaultFs = 50;
defaultN = 5;
defaultall_ages = [7, 14, 21];
defaultbout_rate = [ 1.65,  1.4, 1.4]; % for 7,14 and 21 dpf
defaultSpeed = [ 1.5 , 1.5 , 1.2 ];
defaultBL = [ 0.4 , 0.5 , 0.8 ];
defaultArean_rad = [6.5, 9.2, 12.6];
defaultT = 600; % in sec
defaultOne_side_blindang = 15;  % in degrees s.t. 2*blindangle is the full angle
defaultWall_th = 2;
defaultPLOT = 1;

% parse
vars = inputParser;
addParameter(vars,'Fs',defaultFs);
addParameter(vars,'N',defaultN);
addParameter(vars,'PLOT',defaultPLOT);
addParameter(vars,'all_bout_rate',defaultbout_rate);
addParameter(vars,'all_avg_speed',defaultSpeed);
addParameter(vars,'all_BL',defaultBL);
addParameter(vars,'all_arean_rad',defaultArean_rad);
addParameter(vars,'T',defaultT);
addParameter(vars,'blindang',defaultOne_side_blindang);
addParameter(vars,'wall_th',defaultWall_th);
addParameter(vars,'all_ages',defaultall_ages);


parse(vars,varargin{:})
% make variables
all_bout_rate = vars.Results.all_bout_rate;
all_avg_speed = vars.Results.all_avg_speed;
all_BL = vars.Results.all_BL;
all_arena_rad = vars.Results.all_arean_rad;
T = vars.Results.T;
blindang = vars.Results.blindang;
WD_th = vars.Results.wall_th; %wall distance threshold
PLOT = vars.Results.PLOT;
Fs = vars.Results.Fs;
N = vars.Results.N;
all_ages = vars.Results.all_ages;

% translate time to frame
T = T*Fs; 
dt = 1/Fs'; % delta T

BL = all_BL(age==all_ages); % BL in cm

% fraction from analytical bout rate due to duration of bout exectuation
% [see file AnalyticalVsActualBoutRates.m
est_dec_in_rate = 0.55; %increass to compensate for 'down times' when fush cant bout (found from graphs)
bout_prob = all_bout_rate(age==all_ages)/Fs/est_dec_in_rate; % this is to acount for 'down time';

avg_speed = all_avg_speed(age==all_ages);

% define horizontal and vertical sizes
horz_size = BL;
vert_size = BL/2;

% arena radius
arena_rad = all_arena_rad(age==all_ages);

% REMOVE THIS:
% loadpath = ['Y:\Simulating free swimming larva\'];
% define turning response according to age and type [VR vs group
% experiment]

if strcmp(int_type,'7dpf_exp') % 7 dpf group experiemnt
    %     load([loadpath,'turn_on_ref_diff_WT_AB5fish_fit']);
    %     load([loadpath,'turn_on_ref_diff_Lightallage7fit']);
    load([loadpath,'turn_7dpf_exp']);
    
    turn_fn = @(x) x.*blin+0.5;
    
elseif strcmp(int_type,'14dpf_exp') % 14 dpf group experiemnt
    
    %     load([loadpath,'turn_on_ref_diff_ABallage14_fit']);
    %     load([loadpath,'turn_on_ref_diff_Lightallage14_fit']);
    load([loadpath,'turn_14dpf_exp']);
   
    turn_fn = slm; 
    
elseif strcmp(int_type,'21dpf_exp') % 21 dpf group experiemnt
    %     load([loadpath,'turn_on_ref_diff_ABallage21_fit']);
    %     load([loadpath,'turn_on_ref_diff_Lightallage21_fit']);
    load([loadpath,'turn_21dpf_exp']);
    
    turn_fn = slm; 

elseif strcmp(int_type,'7dpf_VR')
    % set horz and vertical sizes of the fish
    
    
    %     load([loadpath,'mean_turn_prob_by_eye_diff_AB_Age_7_fit_7.mat']);
    load([loadpath,'turn_7dpf_VR']);
    
    p1 = FIT.p1;
    p0 = FIT.p2;
    turn_fn = @(x) x.*p1+p0;
    
elseif strcmp(int_type,'14dpf_VR')
    
    %     load([loadpath,'mean_turn_prob_by_eye_diff_AB_Age_14_fit_8.mat']);
    load([loadpath,'turn_14dpf_VR']);
    turn_fn = slm; % linear on eye diff
    
elseif strcmp(int_type,'21dpf_VR')
    
    load([loadpath,'turn_21dpf_VR']);
    turn_fn = slm; % linear on eye diff
elseif strcmp(int_type,'non_social')
    turn_fn = @(x) x.*0 + 0.5; % no reponse to neighbors always outputs 0.5
    
end

load([loadpath,'wall_response']);

wall_fun = @(x) aexp.*exp(bexp.*x); % note - the function is for BL not cm to match all age groups


%%

% turning parameters:
Wall_sig = 45; % % sigma of normal distribution for turning angles from wall
% WD_th = 2; % BL (since it is age sependend) below this distance fish respond to walls

Neigh_sig = 30; % sigma for normal distribution turning angles from neighbros

% probability to stop chaining turns
Pchange_turn = 1-0.5;

% define a sigmoidal bout shape (for nice visualizaiton of motion)
xsig = (0:7)/Fs;
x0 = 3.5/Fs;
maxV = avg_speed; % taken from data of 7 dpf (we may need to change for the others)
slope = 41; % set the slope (or add a distribution)


% define arena size and bounderies:
rad = arena_rad/2; % arena radius in cm
rad_start = rad; % radius to start agents in
boundx = sind(0:360)*rad;
boundy = cosd(0:360)*rad;


%% run simulaitons:

% define variables:
% define state of each fish (0 - idle, 1 - in bout)
State = zeros(N,T);

% variables for position
x = ones(N,T);
y = ones(N,T);

% random starting positions:
[x_strt,y_strt] = dotsInAcircle(N,0.9*(rad_start),(rad_start)/3,[0 0]);

x(:,1) = x_strt;
y(:,1) = y_strt;
% relPosX = zeros(N,T,N);
% relPosY = zeros(N,T,N);

% variables for neighbor relative angle, orientation and distance
relAng = zeros(N,T,N);
relOri = zeros(N,T,N);
relDist = zeros(N,T,N);

% fish body angle
angle = zeros(N,T);

% startig angle
angle(:,1) = rand(N,1)*360; % set random angles for start (angle is angle from north)

% velocity vectors
Vx = zeros(N,T);
Vy = zeros(N,T);

% speed
Speed = zeros(N,T);

% distance to wall
wallD = zeros(N,T);


% vector from center of arena to fish
Norm = calculateNorm([x(:,1) y(:,1)]);
% distance from wall:
wallD(:,1) = (rad - Norm)/BL; % transform to body length since wall distance funciton in in BL

% save the angle of fish direction from the norm vec from fish to center of
% the wall. values < 180 means wall is to the left (90 is perpendicular),
% values > 180 means wall is to the right (270 is perpendicular)
angFromWall = zeros(N,T);

% normalized vector from center towards fish:
vec_2_fish = [x(:,1)./Norm y(:,1)./Norm];

% starting angle from wall as the angle between the norm to the
% fish and the angle of motion. This gives 0 360 angle. ang < 180
% means wall is to the left, ang > 180 wall is to the right (THIS
% Was double checked)
[~,ang_from_norm] = angOfVectors([sind(angle(:,1)) cosd(angle(:,1))],vec_2_fish);
angFromWall(:,1) = ang_from_norm;

Bouts = zeros(N,T); % var for bouts
all_p = zeros(N,T); % var for proobabilities

% variable for turning left or right, used when no interactions
% exist 1 is , -1 is


% variable for the obtained retina angles:
all_pct_ret = zeros(N,T,2,N-1);

from_bout_count = ones(N,1); % counter for time from previous bout;

% claculate relative propertied for all fish
t = 1;
for f = 1:N
    [relDist(f,t,:),relAng(f,t,:),relOri(f,t,:)] = ....
        relativeNeighborProp(f,x(:,t),y(:,t),angle(:,t));
end

% plot(X(:,1:2)',Y(:,1:2)','.');hold on;
% quiver(X(:,1),Y(:,1),sind(angle(:,1)),cosd(angle(:,1)),0.5);
%
% plot(boundx,boundy,'k');
%         p = 0; %
% rng(1)

for t = 2:T-15 % for all points
    
    idle = find(State(:,t)==0); % fish that are not currently moving
    
    for f = 1:length(idle) % loop over idle fish


        % get wall distance:
        tempWD = wallD(idle(f),t-1);
        
        % estimate retina angle using distance and orientations:
        tempD = squeeze(relDist(idle(f),t-1,:));
        neighD = tempD; % neighbor distnace
        neighD(idle(f)) = [];
        
        tempA = squeeze(relAng(idle(f),t-1,:));
        neighA = tempA; % angle to neighbor
        neighA(idle(f)) = [];
        neighA(neighA > 180) = neighA(neighA>180) - 360;
        
        tempO = squeeze(relOri(idle(f),t-1,:));
        neighO = tempO; % neighbor orientations
        neighO(idle(f)) = [];
        
        
        % calcualte retina angle according to horizontal or vertical size
        if strcmp(int_type,'7dpf_exp') || strcmp(int_type,'14dpf_exp')...
                || strcmp(int_type,'14dpf_exp')
            
            % use only fish that are not on the 
            ii_not_in_blind =  neighA < (180-blindang) & neighA > (-180+blindang);
            left_right_side_fish = ones(size(neighA));
            left_right_side_fish(neighA < 0) = -1; % neighbors to the left -1, on the right 1
            left_side_ii = left_right_side_fish<0;
            
            [ret_angle, right_side] = calcVisualAngle(neighD, neighA, neighO,...
                horz_size);
            
%             ret_angle = ret_angle*5;
%              ret_angle = 2 * atand(horz_size/2./neighD);

            
            % calc ret angles according to horzinotal size, distance and
            % relative angle:
            
            pctr = nansum(ret_angle(right_side==1 & ii_not_in_blind));
            pctl = nansum(ret_angle(right_side==0 & ii_not_in_blind)); % neighbors to the left
            
            % save to var retina occupency:
            all_pct_ret(idle(f),t,1,left_side_ii & ii_not_in_blind) = ...
                ret_angle(left_side_ii & ii_not_in_blind);
            all_pct_ret(idle(f),t,2,~left_side_ii & ii_not_in_blind) = ...
                ret_angle(~left_side_ii & ii_not_in_blind);
            
            % difference between the eyes
            Diff = pctr - pctl;
                       
            % calculate probability to turn right from visual input
            if age== 7 % if this is a linear funciton
                predicted_pright = turn_fn(Diff);
                    
            else
                predicted_pright = slmeval(Diff,slm);
            end
            
            
        else
            % use the simple distance to fish to decide the retinal clutter
            ret_angle = 2 * atand((vert_size/2)./neighD);
            
            % calculate predicted response for each fish
            % according to its size and side
            left_right_side_fish = ones(size(neighA));
            left_right_side_fish(neighA < 0) = -1; % neighbors to the left -1, on the right 1
            left_side_ii = left_right_side_fish<0;
            
            % remove angles of fish in the blind spot:
%             blind_spot_ii = neighA > (180-blindang) | neighA < (-180+blindang);
            ii_not_in_blind =  (neighA < (180-blindang) & neighA > (-180+blindang));
            
            
            % save to var retina occupency:
            all_pct_ret(idle(f),t,1,left_side_ii & ii_not_in_blind) = ...
                ret_angle(left_side_ii & ii_not_in_blind);
            all_pct_ret(idle(f),t,2,~left_side_ii & ii_not_in_blind) = ...
                ret_angle(~left_side_ii & ii_not_in_blind);
           
            
            if age ==7 % age 7 has a linear fit 
                temp_pright = turn_fn(ret_angle.*left_right_side_fish*-1); 
            else
                temp_pright = slmeval(ret_angle,turn_fn);
                temp_pright(~left_side_ii) = 1 - temp_pright(~left_side_ii);
            end
            
            % we take a weighted avg on each side (reponse
            % tendency is the difference from 0.5)
            if age == 7 % use a weighted average
                left_response = sum(temp_pright(left_side_ii & ii_not_in_blind).*...
                    ret_angle(left_side_ii &  ii_not_in_blind))./...
                    sum(ret_angle(left_side_ii & ii_not_in_blind))-0.5;
                
                % we take a weighted avg on each side
                right_response = sum(temp_pright(~left_side_ii & ii_not_in_blind).*...
                    ret_angle(~left_side_ii & ii_not_in_blind))./...
                    sum(ret_angle(~left_side_ii & ii_not_in_blind))-0.5;
            else % use a simple average
                left_response = sum(temp_pright(left_side_ii & ii_not_in_blind))./...
                    sum((left_side_ii & ii_not_in_blind))-0.5;
                
                % we take a weighted avg on each side
                right_response = sum(temp_pright(~left_side_ii & ii_not_in_blind))./...
                    sum((~left_side_ii & ii_not_in_blind))-0.5;
                % ADD A SIMPLE AVERAGE:
            end
            
            % we sum up the signed responses on both sides:
            % (nans are when no neighbors) are on one side
            
            predicted_pright = nansum([left_response,right_response])+0.5;
            all_p(idle(f),t) = predicted_pright;
            
        end
            

        % define prob of bout:
        p_bout = bout_prob; % define basic bout prob
        
        % flip a coin to decide on a bout

        xx(f,t) = rand(1);
        bout_flag = xx(f,t) < p_bout;
        Bouts(idle(f),t) = bout_flag;
        
        % if a bout is starting - choose kinematics
        if bout_flag
            from_bout_count(idle(f)) = 1; % reset the counter
            
            % Kinematics phase:
            % for fish switching to bouts, choose angle, and speed and duration
            % according to neighbors and walls.
            
            % turn away from wall as first priority
            
            if tempWD < WD_th
                
                % probability to resppnd to the wall
                pwall = rand(1) < wall_fun(tempWD)+0.5;
                
                if pwall % if fish is turning from wall
                    
                    % if wall is to the left (<180) turn right) and vice
                    % versa (THIS WAS CHECKED and is correct)
                    p_right = angFromWall(idle(f),t-1) < 180;
                    ang = randn(1)*Wall_sig;
                    
                    % choose turn direction according to p_right:
                    ang = abs(ang)*p_right-abs(ang)*(1-p_right);
                    
                    ang = angle(idle(f),t-1) + ang;
                    
                end
            else
                pwall = 0;
            end
            
            % if fish is not responding to wall it is responding to neighbors:
            if pwall==0
%                 
%                 predicted_pright
%                 t
                % calculate the prob to turn left or right:
                p_right = rand(1) < predicted_pright;
                
                ang = randn(1)*Neigh_sig; % angle from a dist with sigma = 30;
                %
                
                % choose turn direction according to p_right:
                % positive is a right turn negative is a left turn
                ang = abs(ang)*p_right-abs(ang)*(1-p_right);
                
                %                 disp(['fish: ',num2str(idle(f)),' neighbor dep ang: ',...
                %                     num2str(ang,3)]);
                %                 pause();
                % change angle from distribution, according to neighbor
                % properties:
                ang = angle(idle(f),t-1) + ang;
                %                         p = p+1;
          
           
            end
            
            % set kinematics
            
            % speed in body length per second
            Sup = maxV./(1+exp(-slope*(xsig-x0)));
            Sdown = maxV./(1+exp(slope*(xsig-x0)));
            
            % trnasform to cm per frame (so we can sum it)
            S = ([Sup Sdown]*BL)/Fs;
            % velocity
            vx = S*sind(ang);
            vy = S*cosd(ang);
            
            % position
            xx = cumsum(vx);
            yy = cumsum(vy);
            
            % check if the expected trajectory will go outside the arena
            dd = calculateNorm([[xx+x(idle(f),t-1)]' [yy+y(idle(f),t-1)]']); % radius
            
            outside = dd > rad; % position outside
            % if needed correct to stop at the wall
            if sum(outside)>0
                vx(outside) = 0;
                vy(outside) = 0;
                xx = cumsum(vx);
                yy = cumsum(vy);
                S(outside) = 0;
            end
            
            % Update phase:
            % update all positions and kinematic variables of moving fish
            State(idle(f),t:t+length(S)-2) = ~outside(1:end-1); % set bout times
            State(idle(f),t:t+length(S)-1) = 1;
            angle(idle(f),t:t+length(S)-1) = ang;
            
            Speed(idle(f),t:t+length(S)-1) = S;
            Vx(idle(f),t:t+length(vx)-1) = vx;
            Vy(idle(f),t:t+length(vy)-1) = vy;
            
            x(idle(f),t:t+length(xx)-1) = xx+x(idle(f),t-1);
            y(idle(f),t:t+length(xx)-1) = yy+y(idle(f),t-1);
            
            
        else
            x(idle(f),t) = x(idle(f),t-1);
            %                     [x(idle(f),t),x(idle(f),t-1)]
            y(idle(f),t) = y(idle(f),t-1);
            angle(idle(f),t) = angle(idle(f),t-1);
            % increase the counter if no bout was performed
            from_bout_count(idle(f)) = from_bout_count(idle(f))+1;
        end
        
        
        
    end
    
    % calculate relative position and angle for each fish
    
    for f = 1:N
        [relDist(f,t,:),relAng(f,t,:),relOri(f,t,:)] = ....
            relativeNeighborProp(f,x(:,t),y(:,t),angle(:,t));
    end
    
    % wall distance
    wallD(:,t) = (rad - calculateNorm([x(:,t) y(:,t)]))/BL;
    
    % vector from center of arena to fish
    Norm = calculateNorm([x(:,t) y(:,t)]);
    vec_2_fish = [x(:,t)./Norm y(:,t)./Norm];
    
    % set starting positions:
    [~,angFromWall(:,t)] = angOfVectors([sind(angle(:,t)) cosd(angle(:,t))],vec_2_fish);
    
    
    if PLOT && t>20
        %                 plot(x(2,1:500)); hold on;
        %                 plot(x(2,1:t));
        %                 hold off;
        plot(x(:,t-20:t)',y(:,t-20:t)','.'); hold on;
        %                 quiver(x(:,t),y(:,t),Vx(:,t),Vy(:,t),0.1);
        plot(boundx,boundy,'k');
        hold off
        axis image;
        %     axis([-500 500 -500 500]);
        %     [X(:,t-1:t+1) State(:,t-1:t+1)]
        title(round(t/Fs,3));
        pause(0.01);
    end
    
    
    if mod(t,100)==0
        disp(t)
    end
    
    
end




