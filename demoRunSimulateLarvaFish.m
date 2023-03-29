% This is a demo explaining how how to run the 'SimulateLarvaFish.m' code
% %used in Harpaz et al, ‘Precise visuomotor transformations underlying collective behavior in larval zebrafish’. 
% 
% In order to run the code, the user needs to set the following parameters:
% 
% age – age of the simulated fish 7 ,14 or 21. 
% 
% Load_path – exact path where the .mat files containing the response functions of the fish are stored [files can be found at: ]
% 
% type – user should uncomment one of the following
% -	age+_VR – for virtual reality based interactions
% -	age+_exp – for experiment based interactions
% -	Non_social – for no interactions.  
% 
% trial_num – a number identifier for the specific repetition (effect the random number generator, for reproducibility). 
% 
% Optional inputs [default values can be found in the Methods section]:
% Fs - frames per sec 
% N - number of fish 
% bout_rate - in Hz sec 
% arena_rad - in cm.
% T - total simulation time [frames]
% BL - fish body length [in cm]
% PLOT - 1-yes(default), 0-no

% Output: x,y - position [NxT matrices]
%         Vx,Vy - velocitiy [NxT matrices]
%         speed - norm of velocity [NxT]
%         angle - heading angle [in degrees]
%         State - indicates which behavior was axecuted ....
%         WallD - distance to closes wall [in cm];





%...........Local Variable definitions..........
age = 7; % 7, 14 or 21

% 
load_path = ['/path_to_folder_with_dot_mat_files/'];


type = sprintf('%ddpf_VR',age); % for interactions based on virtual assay
% type = sprintf('%ddpf_exp',age); % for interactions based on group experiments

% type = 'non_social'; % for no socaial interactions

trial_num = 4; % simulation repetition
rng(trial_num); % for reproducability

%.................Main Function.................
[x,y,Vx,Vy,Speed, angle, State, wallD] = ...
    SimulateLarvaFishGitHub(type, age, load_path, 'PLOT', 0);








