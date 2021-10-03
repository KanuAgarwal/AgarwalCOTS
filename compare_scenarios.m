%% SETUP
clear all

% MODEL ===================================================================
% DATA --------------------------------------------------------------------
% Get connectivity matrices
load IdentifyKeySources/ConnectivityMatrices_Model_A_2002_P7

% Load in Australian map outline
load AustOutline

% Load in latitude longitude data
load IdentifyKeySources/original_centroids

% Rename variables in lat and long
lat = lg;
lon = lt;
clear lg lt


% PARAMETERS --------------------------------------------------------------
% How long do we want to run the simulation for
t_end = 50;                     % time in years

% Get the number of reefs from the lat, long data
num_reefs = length(lat);

% Parameters struct: store all parameters in a struct
% Constant parameter values for all reefs
% Estimated by Morello et al. (2014)
params.p_tilde = 0.258;         % effect of fast-growing coral on COTS     
params.M_cots = 2.56;           % natural mortality of COTS 
params.p_1_f = 0.129;           % effect of COTS on fast-growing coral

% Known or arbitrarily chosen by Morello et al. (2014)
params.r_f = 0.5;               % intrinsic growth rate of fast-growing coral
params.K_f = 1;                 % carrying capacity of fast-growing coral
params.p_2_f = 10;              % effect of COTS on fast-growing coral

% Known or arbitrarily chosen by me
params.r_c = 0.1;               % coral larvae reproduction rate
params.r_s = 5000;              % starfish larvae reproduction rate

% Connectivity matrices from Bode et al. (2012)
params.omega_c = psurv_d02_1122_P7;             % coral
params.omega_s = psurv_d02_1122_P7;             % starfish


% INITIAL SYSTEM STATE ----------------------------------------------------
% CORAL
% Percentage of fast-growing coral = 80% everywhere
initial_state.C_0_f = 0.8 * params.K_f * ones(num_reefs, 1);

% STARFISH
% Number of COTS aged 2+
initial_state.N_0_2 = zeros(num_reefs, 1);

% Look for reefs within the initiation box, and put some starfish there
for i = 1:num_reefs
    if (lon(i) > -17 && lon(i) < -14.75) && (lat(i) > 145 && lat(i) < 147)
        initial_state.N_0_2(i) = 100;
    end
end

% Initialise age 1 and age 0 COTS based on Morello initial conditions
initial_state.N_0_1 = initial_state.N_0_2 * exp(params.M_cots);
initial_state.N_0_0 = initial_state.N_0_2 * exp(2*params.M_cots);


%% SCENARIO 0: No control, but simulation run with same conditions

% CONTROL EFFORT ----------------------------------------------------------
% No control
control_effort_s0 = 0;


% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s0, C_y_f_s0, N_y_2_s0, N_y_1_s0, N_y_0_s0] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s0);

% Calculate coral cover and cots over time
coral_s0 = sum(C_y_f_s0, 1);
starfish_age2_s0 = sum(N_y_2_s0, 1);
starfish_age1_s0 = sum(N_y_1_s0, 1);
starfish_larvae_s0 = sum(N_y_0_s0, 1);

% Calculate coral cover in initiation box over time
coral_box_s0 = zeros(1, length(coral_s0));
starfish_age2_box_s0 = zeros(1, length(starfish_age2_s0));

% Loop over time
for t = 1:t_end+1
    % Initialise a sum to keep track of the coral cover
    coral_sum = 0;
    starfish_sum = 0;
    % Loop over every reef
    for i = 1:num_reefs
        % If the reef is in the initiation box, add the coral cover to sum
        if (lon(i) > -17 && lon(i) < -14.75) && (lat(i) > 145 && lat(i) < 147)
            coral_sum = coral_sum + C_y_f_s0(i, t);
            starfish_sum = starfish_sum + N_y_2_s0(i, t);
        end
    end
    % Assign the sum to the coral cover matrix
    coral_box_s0(t) = coral_sum;
    starfish_age2_box_s0(t) = starfish_sum;
end


%% SCENARIO 1: Control at initiation box for 50 years

% CONTROL EFFORT ----------------------------------------------------------
% Cull all the starfish only in the initiation box for all 50 years
control_effort_s1 = zeros(num_reefs, t_end);
for i = 1:num_reefs
    if (lon(i) > -17 && lon(i) < -14.75) && (lat(i) > 145 && lat(i) < 147)
        control_effort_s1(i, :) = 1;
    end
end

% Count number of reefs controlled in control scenario 1
num_reefs_control_s1 = nnz(control_effort_s1(:, 1))
budget_s1 = sum(control_effort_s1, 'all')


% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s1, C_y_f_s1, N_y_2_s1, N_y_1_s1, N_y_0_s1] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s1);

% Calculate coral cover and cots over time
coral_s1 = sum(C_y_f_s1, 1);
starfish_age2_s1 = sum(N_y_2_s1, 1);
starfish_age1_s1 = sum(N_y_1_s1, 1);
starfish_larvae_s1 = sum(N_y_0_s1, 1);

% Calculate coral cover in initiation box over time
coral_box_s1 = zeros(1, length(coral_s1));
starfish_age2_box_s1 = zeros(1, length(starfish_age2_s1));

% Loop over time
for t = 1:t_end+1
    % Initialise a sum to keep track of the coral cover
    coral_sum = 0;
    starfish_sum = 0;
    % Loop over every reef
    for i = 1:num_reefs
        % If the reef is in the initiation box, add the coral cover to sum
        if (lon(i) > -17 && lon(i) < -14.75) && (lat(i) > 145 && lat(i) < 147)
            coral_sum = coral_sum + C_y_f_s1(i, t);
            starfish_sum = starfish_sum + N_y_2_s1(i, t);
        end
    end
    % Assign the sum to the coral cover matrix
    coral_box_s1(t) = coral_sum;
    starfish_age2_box_s1(t) = starfish_sum;
end


%% SCENARIO 2: Control at initiation box +10% buffer area for 50 years

% CONTROL EFFORT ----------------------------------------------------------
% Cull all starfish in the initiation box + 10% area buffer for 50 years
control_effort_s2 = zeros(num_reefs, t_end);
for i = 1:num_reefs
    if (lon(i) > -17.11 && lon(i) < -14.64) && (lat(i) > 144.75 && lat(i) < 147.25)
        control_effort_s2(i, :) = 1;
    end
end

% Count number of reefs controlled in control scenario 2
num_reefs_control_s2 = nnz(control_effort_s2(:, 1))
budget_s2 = sum(control_effort_s2, 'all')


% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s2, C_y_f_s2, N_y_2_s2, N_y_1_s2, N_y_0_s2] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s2);

% Calculate coral cover and cots over time
coral_s2 = sum(C_y_f_s2, 1);
starfish_age2_s2 = sum(N_y_2_s2, 1);
starfish_age1_s2 = sum(N_y_1_s2, 1);
starfish_larvae_s2 = sum(N_y_0_s2, 1);

% Calculate coral cover in initiation box over time
coral_box_s2 = zeros(1, length(coral_s2));
starfish_age2_box_s2 = zeros(1, length(starfish_age2_s2));

% Loop over time
for t = 1:t_end+1
    % Initialise a sum to keep track of the coral cover
    coral_sum = 0;
    starfish_sum = 0;
    % Loop over every reef
    for i = 1:num_reefs
        % If the reef is in the initiation box, add the coral cover to sum
        if (lon(i) > -17 && lon(i) < -14.75) && (lat(i) > 145 && lat(i) < 147)
            coral_sum = coral_sum + C_y_f_s2(i, t);
            starfish_sum = starfish_sum + N_y_2_s2(i, t);
        end
    end
    % Assign the sum to the coral cover matrix
    coral_box_s2(t) = coral_sum;
    starfish_age2_box_s2(t) = starfish_sum;
end


%% SCENARIO 3: Control at twice the area but half the effort for 50 years

% CONTROL EFFORT ----------------------------------------------------------
% Cull all starfish in the initiation box + 10% area buffer for 50 years
control_effort_s3 = zeros(num_reefs, t_end);
for i = 1:num_reefs
    if (lon(i) > -18.45 && lon(i) < -14) && (lat(i) > 143 && lat(i) < 147.75)
        control_effort_s3(i, :) = 0.5;
    end
end

% Print the total number of reefs and 'budget' for control scenario
num_reefs_control_s3 = nnz(control_effort_s3(:, 1))
budget_s3 = sum(control_effort_s3, 'all')


% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s3, C_y_f_s3, N_y_2_s3, N_y_1_s3, N_y_0_s3] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s3);

% Calculate coral cover and cots over time
coral_s3 = sum(C_y_f_s3, 1);
starfish_age2_s3 = sum(N_y_2_s3, 1);
starfish_age1_s3 = sum(N_y_1_s3, 1);
starfish_larvae_s3 = sum(N_y_0_s3, 1);

% Calculate coral cover in initiation box over time
coral_box_s3 = zeros(1, length(coral_s3));
starfish_age2_box_s3 = zeros(1, length(starfish_age2_s3));

% Loop over time
for t = 1:t_end+1
    % Initialise a sum to keep track of the coral cover
    coral_sum = 0;
    starfish_sum = 0;
    % Loop over every reef
    for i = 1:num_reefs
        % If the reef is in the initiation box, add the coral cover to sum
        if (lon(i) > -17 && lon(i) < -14.75) && (lat(i) > 145 && lat(i) < 147)
            coral_sum = coral_sum + C_y_f_s3(i, t);
            starfish_sum = starfish_sum + N_y_2_s3(i, t);
        end
    end
    % Assign the sum to the coral cover matrix
    coral_box_s3(t) = coral_sum;
    starfish_age2_box_s3(t) = starfish_sum;
end


%% SCENARIO 4: Control quadruple the area but same total effort for 50 years

% CONTROL EFFORT ----------------------------------------------------------
% Cull all starfish in the initiation box + 10% area buffer for 50 years
control_effort_s4 = zeros(num_reefs, t_end);
for i = 1:num_reefs
    if (lon(i) > -20 && lon(i) < -12.69) && (lat(i) > 142 && lat(i) < 152)
        control_effort_s4(i, :) = 0.25;
    end
end

% Print the total number of reefs and 'budget' for control scenario
num_reefs_control_s4 = nnz(control_effort_s4(:, 1))
budget_s4 = sum(control_effort_s4, 'all')


% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s4, C_y_f_s4, N_y_2_s4, N_y_1_s4, N_y_0_s4] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s4);

% Calculate coral cover and cots over time
coral_s4 = sum(C_y_f_s4, 1);
starfish_age2_s4 = sum(N_y_2_s4, 1);
starfish_age1_s4 = sum(N_y_1_s4, 1);
starfish_larvae_s4 = sum(N_y_0_s4, 1);

% Calculate coral cover in initiation box over time
coral_box_s4 = zeros(1, length(coral_s4));
starfish_age2_box_s4 = zeros(1, length(starfish_age2_s4));

% Loop over time
for t = 1:t_end+1
    % Initialise a sum to keep track of the coral cover
    coral_sum = 0;
    starfish_sum = 0;
    % Loop over every reef
    for i = 1:num_reefs
        % If the reef is in the initiation box, add the coral cover to sum
        if (lon(i) > -17 && lon(i) < -14.75) && (lat(i) > 145 && lat(i) < 147)
            coral_sum = coral_sum + C_y_f_s4(i, t);
            starfish_sum = starfish_sum + N_y_2_s4(i, t);
        end
    end
    % Assign the sum to the coral cover matrix
    coral_box_s4(t) = coral_sum;
    starfish_age2_box_s4(t) = starfish_sum;
end


%% PLOT: Comparison of coral cover over time

% Total coral cover
figure(11), clf, hold on, grid on
plot(t_vec_s0, coral_s0, 'Linewidth', 2)
plot(t_vec_s1, coral_s1, 'Linewidth', 2)
plot(t_vec_s2, coral_s2, 'Linewidth', 2)
plot(t_vec_s3, coral_s3, 'Linewidth', 2)
plot(t_vec_s4, coral_s4, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total coral cover')
title('Total coral cover on GBR for different control scenarios')
legend('No control', '168 reefs, 8400 budget', '185 reefs, 9250 budget', ...
       '336 reefs, 8400 budget', '672 reefs, 8400 budget', ...
       'Location', 'NorthEastOutside')

% Coral cover at initiation box
figure(12), clf, hold on, grid on
plot(t_vec_s0, coral_box_s0, 'Linewidth', 2)
plot(t_vec_s1, coral_box_s1, 'Linewidth', 2)
plot(t_vec_s2, coral_box_s2, 'Linewidth', 2)
plot(t_vec_s3, coral_box_s3, 'Linewidth', 2)
plot(t_vec_s4, coral_box_s4, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total coral cover')
title('Total coral cover at initiation box for different control scenarios')
legend('No control', '168 reefs, 8400 budget', '185 reefs, 9250 budget', ...
       '336 reefs, 8400 budget', '672 reefs, 8400 budget', ...
       'Location', 'NorthEastOutside')

   
%% PLOT: Comparison of starfish population over time

% GBR
figure(13), clf, hold on, grid on
plot(t_vec_s0, starfish_age2_s0, 'Linewidth', 2)
plot(t_vec_s1, starfish_age2_s1, 'Linewidth', 2)
plot(t_vec_s2, starfish_age2_s2, 'Linewidth', 2)
plot(t_vec_s3, starfish_age2_s3, 'Linewidth', 2)
plot(t_vec_s4, starfish_age2_s4, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total age 2+ starfish population')
title('Total starfish population on GBR for different control scenarios')
legend('No control', '168 reefs, 8400 budget', '185 reefs, 9250 budget', ...
       '336 reefs, 8400 budget', '672 reefs, 8400 budget', ...
       'Location', 'NorthEastOutside')
   

% GBR - control only
figure(14), clf, hold on, grid on
plot(t_vec_s1, starfish_age2_s1, 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980])
plot(t_vec_s2, starfish_age2_s2, 'Linewidth', 2, 'Color', [0.9290 0.6940 0.1250])
plot(t_vec_s3, starfish_age2_s3, 'Linewidth', 2, 'Color', [0.4940 0.1840 0.5560])
plot(t_vec_s4, starfish_age2_s4, 'Linewidth', 2, 'Color', [0.4660 0.6740 0.1880])
xlabel('Time (years)')
ylabel('Total age 2+ starfish population')
title('Total starfish population on GBR for different control scenarios')
legend('168 reefs, 8400 budget', '185 reefs, 9250 budget', ...
       '336 reefs, 8400 budget', '672 reefs, 8400 budget', ...
       'Location', 'NorthEastOutside')


% Initiation box
figure(15), clf, hold on, grid on
plot(t_vec_s0, starfish_age2_box_s0, 'Linewidth', 2)
plot(t_vec_s1, starfish_age2_box_s1, 'Linewidth', 2)
plot(t_vec_s2, starfish_age2_box_s2, 'Linewidth', 2)
plot(t_vec_s3, starfish_age2_box_s3, 'Linewidth', 2)
plot(t_vec_s4, starfish_age2_box_s4, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total age 2+ starfish population')
title('Total starfish population at initiation box for different control scenarios')
legend('No control', '168 reefs, 8400 budget', '185 reefs, 9250 budget', ...
       '336 reefs, 8400 budget', '672 reefs, 8400 budget', ...
       'Location', 'NorthEastOutside')