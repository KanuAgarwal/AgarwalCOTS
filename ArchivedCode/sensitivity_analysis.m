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
lon = lg;
lat = lt;
clear lg lt


% PARAMETERS --------------------------------------------------------------
% How long do we want to run the simulation for
t_end = 100;                     % time in years

% Get the number of reefs from the lat, long data
num_reefs = length(lon);

% Parameters struct: store all parameters in a struct
% Constant parameter values for all reefs
% Estimated by Morello et al. (2014)
params.p_tilde = 0.258;         % effect of fast-growing coral on COTS     
params.M_cots = 2.56;           % natural mortality of COTS 
params.p_1_f = 0.129/2500;      % effect of COTS on fast-growing coral

% Known or arbitrarily chosen by Morello et al. (2014)
params.r_f = 0.5;               % intrinsic growth rate of fast-growing coral
params.K_f = 1;                 % carrying capacity of fast-growing coral
params.p_2_f = 10/2500;         % effect of COTS on fast-growing coral

% Known or arbitrarily chosen by me
params.r_c = 0.1;               % coral larvae reproduction rate
params.r_s = 5000;              % starfish larvae reproduction rate

% Connectivity matrices from Bode et al. (2012)
params.omega_c = psurv_d02_1122_P7;             % coral
params.omega_s = psurv_d02_1122_P7;             % starfish

% Latitude and longitude for starfish larval calculation
params.lon = lat;
params.lat = lon;

% Using metapopulation model equation for larval dispersal
dispersal_eq = 1;


%% CONTROL SCENARIOS

% SCENARIO 0: No control --------------------------------------------------
control_effort_s0 = 0;

% SCENARIO 1: Control at initiation box -----------------------------------
control_effort_s1 = zeros(num_reefs, t_end);
for i = 1:num_reefs
    if (lat(i) > -17 && lat(i) < -14.75) && (lon(i) > 145 && lon(i) < 147)
        control_effort_s1(i, :) = 1;
    end
end

% Count number of reefs controlled in control scenario 1
num_reefs_control_s1 = nnz(control_effort_s1(:, 1));
budget_s1 = sum(control_effort_s1, 'all');

% SCENARIO 2: Control at twice the area but half the effort ---------------
control_effort_s2 = zeros(num_reefs, t_end);
for i = 1:num_reefs
    if (lat(i) > -18.45 && lat(i) < -14) && (lon(i) > 143 && lon(i) < 147.75)
        control_effort_s2(i, :) = 0.5;
    end
end

% SCENARIO 3: Control four times the area but a quarter of the effort -----
control_effort_s3 = zeros(num_reefs, t_end);
for i = 1:num_reefs
    if (lat(i) > -20 && lat(i) < -12.69) && (lon(i) > 142 && lon(i) < 152)
        control_effort_s3(i, :) = 0.25;
    end
end

% SCENARIO 4: Control entire reef based on total budget -------------------
effort_per_year = budget_s1 / t_end;
effort_per_reef = effort_per_year / num_reefs;
control_effort_s4 = effort_per_reef * ones(num_reefs, t_end);


%% CASE 1: 5 COTS AT EACH REEF ============================================
% INITIAL SYSTEM STATE ----------------------------------------------------
% CORAL
% Percentage of fast-growing coral = 80% everywhere
initial_state.C_0_f = 0.8 * params.K_f * ones(num_reefs, 1);

% STARFISH
% Number of COTS aged 2+
initial_state.N_0_2 = zeros(num_reefs, 1);

% Look for reefs within the initiation box, and put some starfish there
for i = 1:num_reefs
    if (lat(i) > -17 && lat(i) < -14.75) && (lon(i) > 145 && lon(i) < 147)
        initial_state.N_0_2(i) = 5;
    end
end

% Initialise age 1 and age 0 COTS based on Morello initial conditions
initial_state.N_0_1 = initial_state.N_0_2 * exp(params.M_cots);
initial_state.N_0_0 = initial_state.N_0_2 * exp(2*params.M_cots);


% CONTROL SCENARIO 0 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c1_s0, C_y_f_c1_s0, N_y_2_c1_s0, N_y_1_c1_s0, N_y_0_c1_s0, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s0, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c1_s0 = sum(C_y_f_c1_s0, 1);
starfish_age2_c1_s0 = sum(N_y_2_c1_s0, 1);
starfish_age1_c1_s0 = sum(N_y_1_c1_s0, 1);
starfish_age0_c1_s0 = sum(N_y_0_c1_s0, 1);

% Calculate coral cover in initiation box over time
[coral_box_c1_s0, starfish_age2_box_c1_s0, starfish_age1_box_c1_s0, starfish_age0_box_c1_s0] = ...
    calculate_population_box(t_end, C_y_f_c1_s0, N_y_2_c1_s0, N_y_1_c1_s0, N_y_0_c1_s0, num_reefs, lon, lat);


% CONTROL SCENARIO 1 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c2_s1, C_y_f_c1_s1, N_y_2_c1_s1, N_y_1_c1_s1, N_y_0_c1_s1, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s1, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c1_s1 = sum(C_y_f_c1_s1, 1);
starfish_age2_c1_s1 = sum(N_y_2_c1_s1, 1);
starfish_age1_c1_s1 = sum(N_y_1_c1_s1, 1);
starfish_age0_c1_s1 = sum(N_y_0_c1_s1, 1);

% Calculate coral cover in initiation box over time
[coral_box_c1_s1, starfish_age2_box_c1_s1, starfish_age1_box_c1_s1, starfish_age0_box_c1_s1] = ...
    calculate_population_box(t_end, C_y_f_c1_s1, N_y_2_c1_s1, N_y_1_c1_s1, N_y_0_c1_s1, num_reefs, lon, lat);


% CONTROL SCENARIO 2 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c1_s2, C_y_f_c1_s2, N_y_2_c1_s2, N_y_1_c1_s2, N_y_0_c1_s2, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s2, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c1_s2 = sum(C_y_f_c1_s2, 1);
starfish_age2_c1_s2 = sum(N_y_2_c1_s2, 1);
starfish_age1_c1_s2 = sum(N_y_1_c1_s2, 1);
starfish_age0_c1_s2 = sum(N_y_0_c1_s2, 1);

% Calculate coral cover in initiation box over time
[coral_box_c1_s2, starfish_age2_box_c1_s2, starfish_age1_box_c1_s2, starfish_age0_box_c1_s2] = ...
    calculate_population_box(t_end, C_y_f_c1_s2, N_y_2_c1_s2, N_y_1_c1_s2, N_y_0_c1_s2, num_reefs, lon, lat);


% CONTROL SCENARIO 3 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c1_s3, C_y_f_c1_s3, N_y_2_c1_s3, N_y_1_c1_s3, N_y_0_c1_s3, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s3, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c1_s3 = sum(C_y_f_c1_s3, 1);
starfish_age2_c1_s3 = sum(N_y_2_c1_s3, 1);
starfish_age1_c1_s3 = sum(N_y_1_c1_s3, 1);
starfish_age0_c1_s3 = sum(N_y_0_c1_s3, 1);

% Calculate coral cover in initiation box over time
[coral_box_c1_s3, starfish_age2_box_c1_s3, starfish_age1_box_c1_s3, starfish_age0_box_c1_s3] = ...
    calculate_population_box(t_end, C_y_f_c1_s3, N_y_2_c1_s3, N_y_1_c1_s3, N_y_0_c1_s3, num_reefs, lon, lat);


% CONTROL SCENARIO 4 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c1_s4, C_y_f_c1_s4, N_y_2_c1_s4, N_y_1_c1_s4, N_y_0_c1_s4, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s4, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c1_s4 = sum(C_y_f_c1_s4, 1);
starfish_age2_c1_s4 = sum(N_y_2_c1_s4, 1);
starfish_age1_c1_s4 = sum(N_y_1_c1_s4, 1);
starfish_age0_c1_s4 = sum(N_y_0_c1_s4, 1);

% Calculate coral cover in initiation box over time
[coral_box_c1_s4, starfish_age2_box_c1_s4, starfish_age1_box_c1_s4, starfish_age0_box_c1_s4] = ...
    calculate_population_box(t_end, C_y_f_c1_s4, N_y_2_c1_s4, N_y_1_c1_s4, N_y_0_c1_s4, num_reefs, lon, lat);


%% CASE 2: 10 COTS AT EACH REEF ===========================================
% INITIAL SYSTEM STATE ----------------------------------------------------
% CORAL
% Percentage of fast-growing coral = 80% everywhere
initial_state.C_0_f = 0.8 * params.K_f * ones(num_reefs, 1);

% STARFISH
% Number of COTS aged 2+
initial_state.N_0_2 = zeros(num_reefs, 1);

% Look for reefs within the initiation box, and put some starfish there
for i = 1:num_reefs
    if (lat(i) > -17 && lat(i) < -14.75) && (lon(i) > 145 && lon(i) < 147)
        initial_state.N_0_2(i) = 10;
    end
end

% Initialise age 1 and age 0 COTS based on Morello initial conditions
initial_state.N_0_1 = initial_state.N_0_2 * exp(params.M_cots);
initial_state.N_0_0 = initial_state.N_0_2 * exp(2*params.M_cots);


% CONTROL SCENARIO 0 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c2_s0, C_y_f_c2_s0, N_y_2_c2_s0, N_y_1_c2_s0, N_y_0_c2_s0, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s0, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c2_s0 = sum(C_y_f_c2_s0, 1);
starfish_age2_c2_s0 = sum(N_y_2_c2_s0, 1);
starfish_age1_c2_s0 = sum(N_y_1_c2_s0, 1);
starfish_age0_c2_s0 = sum(N_y_0_c2_s0, 1);

% Calculate coral cover in initiation box over time
[coral_box_c2_s0, starfish_age2_box_c2_s0, starfish_age1_box_c2_s0, starfish_age0_box_c2_s0] = ...
    calculate_population_box(t_end, C_y_f_c2_s0, N_y_2_c2_s0, N_y_1_c2_s0, N_y_0_c2_s0, num_reefs, lon, lat);


% CONTROL SCENARIO 1 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c2_s1, C_y_f_c2_s1, N_y_2_c2_s1, N_y_1_c2_s1, N_y_0_c2_s1, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s1, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c2_s1 = sum(C_y_f_c2_s1, 1);
starfish_age2_c2_s1 = sum(N_y_2_c2_s1, 1);
starfish_age1_c2_s1 = sum(N_y_1_c2_s1, 1);
starfish_age0_c2_s1 = sum(N_y_0_c2_s1, 1);

% Calculate coral cover in initiation box over time
[coral_box_c2_s1, starfish_age2_box_c2_s1, starfish_age1_box_c2_s1, starfish_age0_box_c2_s1] = ...
    calculate_population_box(t_end, C_y_f_c2_s1, N_y_2_c2_s1, N_y_1_c2_s1, N_y_0_c2_s1, num_reefs, lon, lat);


% CONTROL SCENARIO 2 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c2_s2, C_y_f_c2_s2, N_y_2_c2_s2, N_y_1_c2_s2, N_y_0_c2_s2, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s2, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c2_s2 = sum(C_y_f_c2_s2, 1);
starfish_age2_c2_s2 = sum(N_y_2_c2_s2, 1);
starfish_age1_c2_s2 = sum(N_y_1_c2_s2, 1);
starfish_age0_c2_s2 = sum(N_y_0_c2_s2, 1);

% Calculate coral cover in initiation box over time
[coral_box_c2_s2, starfish_age2_box_c2_s2, starfish_age1_box_c2_s2, starfish_age0_box_c2_s2] = ...
    calculate_population_box(t_end, C_y_f_c2_s2, N_y_2_c2_s2, N_y_1_c2_s2, N_y_0_c2_s2, num_reefs, lon, lat);


% CONTROL SCENARIO 3 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c2_s3, C_y_f_c2_s3, N_y_2_c2_s3, N_y_1_c2_s3, N_y_0_c2_s3, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s3, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c2_s3 = sum(C_y_f_c2_s3, 1);
starfish_age2_c2_s3 = sum(N_y_2_c2_s3, 1);
starfish_age1_c2_s3 = sum(N_y_1_c2_s3, 1);
starfish_age0_c2_s3 = sum(N_y_0_c2_s3, 1);

% Calculate coral cover in initiation box over time
[coral_box_c2_s3, starfish_age2_box_c2_s3, starfish_age1_box_c2_s3, starfish_age0_box_c2_s3] = ...
    calculate_population_box(t_end, C_y_f_c2_s3, N_y_2_c2_s3, N_y_1_c2_s3, N_y_0_c2_s3, num_reefs, lon, lat);


% CONTROL SCENARIO 4 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c2_s4, C_y_f_c2_s4, N_y_2_c2_s4, N_y_1_c2_s4, N_y_0_c2_s4, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s4, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c2_s4 = sum(C_y_f_c2_s4, 1);
starfish_age2_c2_s4 = sum(N_y_2_c2_s4, 1);
starfish_age1_c2_s4 = sum(N_y_1_c2_s4, 1);
starfish_age0_c2_s4 = sum(N_y_0_c2_s4, 1);

% Calculate coral cover in initiation box over time
[coral_box_c2_s4, starfish_age2_box_c2_s4, starfish_age1_box_c2_s4, starfish_age0_box_c2_s4] = ...
    calculate_population_box(t_end, C_y_f_c2_s4, N_y_2_c2_s4, N_y_1_c2_s4, N_y_0_c2_s4, num_reefs, lon, lat);


%% CASE 3: 25 COTS AT EACH REEF ===========================================
% INITIAL SYSTEM STATE ----------------------------------------------------
% CORAL
% Percentage of fast-growing coral = 80% everywhere
initial_state.C_0_f = 0.8 * params.K_f * ones(num_reefs, 1);

% STARFISH
% Number of COTS aged 2+
initial_state.N_0_2 = zeros(num_reefs, 1);

% Look for reefs within the initiation box, and put some starfish there
for i = 1:num_reefs
    if (lat(i) > -17 && lat(i) < -14.75) && (lon(i) > 145 && lon(i) < 147)
        initial_state.N_0_2(i) = 25;
    end
end

% Initialise age 1 and age 0 COTS based on Morello initial conditions
initial_state.N_0_1 = initial_state.N_0_2 * exp(params.M_cots);
initial_state.N_0_0 = initial_state.N_0_2 * exp(2*params.M_cots);


% CONTROL SCENARIO 0 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c3_s0, C_y_f_c3_s0, N_y_2_c3_s0, N_y_1_c3_s0, N_y_0_c3_s0, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s0, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c3_s0 = sum(C_y_f_c3_s0, 1);
starfish_age2_c3_s0 = sum(N_y_2_c3_s0, 1);
starfish_age1_c3_s0 = sum(N_y_1_c3_s0, 1);
starfish_age0_c3_s0 = sum(N_y_0_c3_s0, 1);

% Calculate coral cover in initiation box over time
[coral_box_c3_s0, starfish_age2_box_c3_s0, starfish_age1_box_c3_s0, starfish_age0_box_c3_s0] = ...
    calculate_population_box(t_end, C_y_f_c3_s0, N_y_2_c3_s0, N_y_1_c3_s0, N_y_0_c3_s0, num_reefs, lon, lat);


% CONTROL SCENARIO 1 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c3_s1, C_y_f_c3_s1, N_y_2_c3_s1, N_y_1_c3_s1, N_y_0_c3_s1, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s1, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c3_s1 = sum(C_y_f_c3_s1, 1);
starfish_age2_c3_s1 = sum(N_y_2_c3_s1, 1);
starfish_age1_c3_s1 = sum(N_y_1_c3_s1, 1);
starfish_age0_c3_s1 = sum(N_y_0_c3_s1, 1);

% Calculate coral cover in initiation box over time
[coral_box_c3_s1, starfish_age2_box_c3_s1, starfish_age1_box_c3_s1, starfish_age0_box_c3_s1] = ...
    calculate_population_box(t_end, C_y_f_c3_s1, N_y_2_c3_s1, N_y_1_c3_s1, N_y_0_c3_s1, num_reefs, lon, lat);


% CONTROL SCENARIO 2 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c3_s2, C_y_f_c3_s2, N_y_2_c3_s2, N_y_1_c3_s2, N_y_0_c3_s2, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s2, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c3_s2 = sum(C_y_f_c3_s2, 1);
starfish_age2_c3_s2 = sum(N_y_2_c3_s2, 1);
starfish_age1_c3_s2 = sum(N_y_1_c3_s2, 1);
starfish_age0_c3_s2 = sum(N_y_0_c3_s2, 1);

% Calculate coral cover in initiation box over time
[coral_box_c3_s2, starfish_age2_box_c3_s2, starfish_age1_box_c3_s2, starfish_age0_box_c3_s2] = ...
    calculate_population_box(t_end, C_y_f_c3_s2, N_y_2_c3_s2, N_y_1_c3_s2, N_y_0_c3_s2, num_reefs, lon, lat);


% CONTROL SCENARIO 3 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c3_s3, C_y_f_c3_s3, N_y_2_c3_s3, N_y_1_c3_s3, N_y_0_c3_s3, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s3, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c3_s3 = sum(C_y_f_c3_s3, 1);
starfish_age2_c3_s3 = sum(N_y_2_c3_s3, 1);
starfish_age1_c3_s3 = sum(N_y_1_c3_s3, 1);
starfish_age0_c3_s3 = sum(N_y_0_c3_s3, 1);

% Calculate coral cover in initiation box over time
[coral_box_c3_s3, starfish_age2_box_c3_s3, starfish_age1_box_c3_s3, starfish_age0_box_c3_s3] = ...
    calculate_population_box(t_end, C_y_f_c3_s3, N_y_2_c3_s3, N_y_1_c3_s3, N_y_0_c3_s3, num_reefs, lon, lat);


% CONTROL SCENARIO 4 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c3_s4, C_y_f_c3_s4, N_y_2_c3_s4, N_y_1_c3_s4, N_y_0_c3_s4, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s4, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c3_s4 = sum(C_y_f_c3_s4, 1);
starfish_age2_c3_s4 = sum(N_y_2_c3_s4, 1);
starfish_age1_c3_s4 = sum(N_y_1_c3_s4, 1);
starfish_age0_c3_s4 = sum(N_y_0_c3_s4, 1);

% Calculate coral cover in initiation box over time
[coral_box_c3_s4, starfish_age2_box_c3_s4, starfish_age1_box_c3_s4, starfish_age0_box_c3_s4] = ...
    calculate_population_box(t_end, C_y_f_c3_s4, N_y_2_c3_s4, N_y_1_c3_s4, N_y_0_c3_s4, num_reefs, lon, lat);


%% CASE 4: 50 COTS AT EACH REEF ===========================================
% INITIAL SYSTEM STATE ----------------------------------------------------
% CORAL
% Percentage of fast-growing coral = 80% everywhere
initial_state.C_0_f = 0.8 * params.K_f * ones(num_reefs, 1);

% STARFISH
% Number of COTS aged 2+
initial_state.N_0_2 = zeros(num_reefs, 1);

% Look for reefs within the initiation box, and put some starfish there
for i = 1:num_reefs
    if (lat(i) > -17 && lat(i) < -14.75) && (lon(i) > 145 && lon(i) < 147)
        initial_state.N_0_2(i) = 50;
    end
end

% Initialise age 1 and age 0 COTS based on Morello initial conditions
initial_state.N_0_1 = initial_state.N_0_2 * exp(params.M_cots);
initial_state.N_0_0 = initial_state.N_0_2 * exp(2*params.M_cots);


% CONTROL SCENARIO 0 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c4_s0, C_y_f_c4_s0, N_y_2_c4_s0, N_y_1_c4_s0, N_y_0_c4_s0, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s0, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c4_s0 = sum(C_y_f_c4_s0, 1);
starfish_age2_c4_s0 = sum(N_y_2_c4_s0, 1);
starfish_age1_c4_s0 = sum(N_y_1_c4_s0, 1);
starfish_age0_c4_s0 = sum(N_y_0_c4_s0, 1);

% Calculate coral cover in initiation box over time
[coral_box_c4_s0, starfish_age2_box_c4_s0, starfish_age1_box_c4_s0, starfish_age0_box_c4_s0] = ...
    calculate_population_box(t_end, C_y_f_c4_s0, N_y_2_c4_s0, N_y_1_c4_s0, N_y_0_c4_s0, num_reefs, lon, lat);


% CONTROL SCENARIO 1 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c4_s1, C_y_f_c4_s1, N_y_2_c4_s1, N_y_1_c4_s1, N_y_0_c4_s1, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s1, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c4_s1 = sum(C_y_f_c4_s1, 1);
starfish_age2_c4_s1 = sum(N_y_2_c4_s1, 1);
starfish_age1_c4_s1 = sum(N_y_1_c4_s1, 1);
starfish_age0_c4_s1 = sum(N_y_0_c4_s1, 1);

% Calculate coral cover in initiation box over time
[coral_box_c4_s1, starfish_age2_box_c4_s1, starfish_age1_box_c4_s1, starfish_age0_box_c4_s1] = ...
    calculate_population_box(t_end, C_y_f_c4_s1, N_y_2_c4_s1, N_y_1_c4_s1, N_y_0_c4_s1, num_reefs, lon, lat);


% CONTROL SCENARIO 2 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c4_s2, C_y_f_c4_s2, N_y_2_c4_s2, N_y_1_c4_s2, N_y_0_c4_s2, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s2, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c4_s2 = sum(C_y_f_c4_s2, 1);
starfish_age2_c4_s2 = sum(N_y_2_c4_s2, 1);
starfish_age1_c4_s2 = sum(N_y_1_c4_s2, 1);
starfish_age0_c4_s2 = sum(N_y_0_c4_s2, 1);

% Calculate coral cover in initiation box over time
[coral_box_c4_s2, starfish_age2_box_c4_s2, starfish_age1_box_c4_s2, starfish_age0_box_c4_s2] = ...
    calculate_population_box(t_end, C_y_f_c4_s2, N_y_2_c4_s2, N_y_1_c4_s2, N_y_0_c4_s2, num_reefs, lon, lat);


% CONTROL SCENARIO 3 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c4_s3, C_y_f_c4_s3, N_y_2_c4_s3, N_y_1_c4_s3, N_y_0_c4_s3, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s3, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c4_s3 = sum(C_y_f_c4_s3, 1);
starfish_age2_c4_s3 = sum(N_y_2_c4_s3, 1);
starfish_age1_c4_s3 = sum(N_y_1_c4_s3, 1);
starfish_age0_c4_s3 = sum(N_y_0_c4_s3, 1);

% Calculate coral cover in initiation box over time
[coral_box_c4_s3, starfish_age2_box_c4_s3, starfish_age1_box_c4_s3, starfish_age0_box_c4_s3] = ...
    calculate_population_box(t_end, C_y_f_c4_s3, N_y_2_c4_s3, N_y_1_c4_s3, N_y_0_c4_s3, num_reefs, lon, lat);


% CONTROL SCENARIO 4 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c4_s4, C_y_f_c4_s4, N_y_2_c4_s4, N_y_1_c4_s4, N_y_0_c4_s4, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s4, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c4_s4 = sum(C_y_f_c4_s4, 1);
starfish_age2_c4_s4 = sum(N_y_2_c4_s4, 1);
starfish_age1_c4_s4 = sum(N_y_1_c4_s4, 1);
starfish_age0_c4_s4 = sum(N_y_0_c4_s4, 1);

% Calculate coral cover in initiation box over time
[coral_box_c4_s4, starfish_age2_box_c4_s4, starfish_age1_box_c4_s4, starfish_age0_box_c4_s4] = ...
    calculate_population_box(t_end, C_y_f_c4_s4, N_y_2_c4_s4, N_y_1_c4_s4, N_y_0_c4_s4, num_reefs, lon, lat);


%% CASE 5: 100 COTS AT EACH REEF ==========================================
% INITIAL SYSTEM STATE ----------------------------------------------------
% CORAL
% Percentage of fast-growing coral = 80% everywhere
initial_state.C_0_f = 0.8 * params.K_f * ones(num_reefs, 1);

% STARFISH
% Number of COTS aged 2+
initial_state.N_0_2 = zeros(num_reefs, 1);

% Look for reefs within the initiation box, and put some starfish there
for i = 1:num_reefs
    if (lat(i) > -17 && lat(i) < -14.75) && (lon(i) > 145 && lon(i) < 147)
        initial_state.N_0_2(i) = 100;
    end
end

% Initialise age 1 and age 0 COTS based on Morello initial conditions
initial_state.N_0_1 = initial_state.N_0_2 * exp(params.M_cots);
initial_state.N_0_0 = initial_state.N_0_2 * exp(2*params.M_cots);


% CONTROL SCENARIO 0 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c5_s0, C_y_f_c5_s0, N_y_2_c5_s0, N_y_1_c5_s0, N_y_0_c5_s0, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s0, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c5_s0 = sum(C_y_f_c5_s0, 1);
starfish_age2_c5_s0 = sum(N_y_2_c5_s0, 1);
starfish_age1_c5_s0 = sum(N_y_1_c5_s0, 1);
starfish_age0_c5_s0 = sum(N_y_0_c5_s0, 1);

% Calculate coral cover in initiation box over time
[coral_box_c5_s0, starfish_age2_box_c5_s0, starfish_age1_box_c5_s0, starfish_age0_box_c5_s0] = ...
    calculate_population_box(t_end, C_y_f_c5_s0, N_y_2_c5_s0, N_y_1_c5_s0, N_y_0_c5_s0, num_reefs, lon, lat);


% CONTROL SCENARIO 1 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c5_s1, C_y_f_c5_s1, N_y_2_c5_s1, N_y_1_c5_s1, N_y_0_c5_s1, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s1, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c5_s1 = sum(C_y_f_c5_s1, 1);
starfish_age2_c5_s1 = sum(N_y_2_c5_s1, 1);
starfish_age1_c5_s1 = sum(N_y_1_c5_s1, 1);
starfish_age0_c5_s1 = sum(N_y_0_c5_s1, 1);

% Calculate coral cover in initiation box over time
[coral_box_c5_s1, starfish_age2_box_c5_s1, starfish_age1_box_c5_s1, starfish_age0_box_c5_s1] = ...
    calculate_population_box(t_end, C_y_f_c5_s1, N_y_2_c5_s1, N_y_1_c5_s1, N_y_0_c5_s1, num_reefs, lon, lat);


% CONTROL SCENARIO 2 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c5_s2, C_y_f_c5_s2, N_y_2_c5_s2, N_y_1_c5_s2, N_y_0_c5_s2, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s2, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c5_s2 = sum(C_y_f_c5_s2, 1);
starfish_age2_c5_s2 = sum(N_y_2_c5_s2, 1);
starfish_age1_c5_s2 = sum(N_y_1_c5_s2, 1);
starfish_age0_c5_s2 = sum(N_y_0_c5_s2, 1);

% Calculate coral cover in initiation box over time
[coral_box_c5_s2, starfish_age2_box_c5_s2, starfish_age1_box_c5_s2, starfish_age0_box_c5_s2] = ...
    calculate_population_box(t_end, C_y_f_c5_s2, N_y_2_c5_s2, N_y_1_c5_s2, N_y_0_c5_s2, num_reefs, lon, lat);


% CONTROL SCENARIO 3 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c5_s3, C_y_f_c5_s3, N_y_2_c5_s3, N_y_1_c5_s3, N_y_0_c5_s3, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s3, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c5_s3 = sum(C_y_f_c5_s3, 1);
starfish_age2_c5_s3 = sum(N_y_2_c5_s3, 1);
starfish_age1_c5_s3 = sum(N_y_1_c5_s3, 1);
starfish_age0_c5_s3 = sum(N_y_0_c5_s3, 1);

% Calculate coral cover in initiation box over time
[coral_box_c5_s3, starfish_age2_box_c5_s3, starfish_age1_box_c5_s3, starfish_age0_box_c5_s3] = ...
    calculate_population_box(t_end, C_y_f_c5_s3, N_y_2_c5_s3, N_y_1_c5_s3, N_y_0_c5_s3, num_reefs, lon, lat);


% CONTROL SCENARIO 4 ------------------------------------------------------
% SOLVE 
% Solve using function which runs simulations
[t_vec_c5_s4, C_y_f_c5_s4, N_y_2_c5_s4, N_y_1_c5_s4, N_y_0_c5_s4, ~] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s4, dispersal_eq);

% CALCULATIONS 
% Calculate coral cover and cots over time
coral_c5_s4 = sum(C_y_f_c5_s4, 1);
starfish_age2_c5_s4 = sum(N_y_2_c5_s4, 1);
starfish_age1_c5_s4 = sum(N_y_1_c5_s4, 1);
starfish_age0_c5_s4 = sum(N_y_0_c5_s4, 1);

% Calculate coral cover in initiation box over time
[coral_box_c5_s4, starfish_age2_box_c5_s4, starfish_age1_box_c5_s4, starfish_age0_box_c5_s4] = ...
    calculate_population_box(t_end, C_y_f_c5_s4, N_y_2_c5_s4, N_y_1_c5_s4, N_y_0_c5_s4, num_reefs, lon, lat);


%% PLOT

% Fontsizes for plotting
axis_FS = 14;
title_FS = 16;
legend_FS = 13;
ticks_FS = 12;

% Colours for plotting
colour_scheme = cbrewer('div', 'RdBu', 4);

% Collate all the data for plotting
initial_starfish_vals = categorical({'5', '10', '25', '50', '100'});
initial_starfish_vals = reordercats(initial_starfish_vals, {'5', '10', '25', '50', '100'});
scenario_0_coral = [coral_c1_s0(end) coral_c2_s0(end) coral_c3_s0(end) coral_c4_s0(end) coral_c5_s0(end)];
scenario_1_coral = [coral_c1_s1(end) coral_c2_s1(end) coral_c3_s1(end) coral_c4_s1(end) coral_c5_s1(end)];
scenario_2_coral = [coral_c1_s2(end) coral_c2_s2(end) coral_c3_s2(end) coral_c4_s2(end) coral_c5_s2(end)];
scenario_3_coral = [coral_c1_s3(end) coral_c2_s3(end) coral_c3_s3(end) coral_c4_s3(end) coral_c5_s3(end)];
scenario_4_coral = [coral_c1_s4(end) coral_c2_s4(end) coral_c3_s4(end) coral_c4_s4(end) coral_c5_s4(end)];

% Plot final coral cover for all scenarios against number of initial adult starfish
figure(1), clf, hold on
scatter(initial_starfish_vals, scenario_0_coral, 100, 'black', 'filled')
scatter(initial_starfish_vals, scenario_1_coral, 100, colour_scheme(1, :), 'filled')
scatter(initial_starfish_vals, scenario_2_coral, 100, colour_scheme(2, :), 'filled')
scatter(initial_starfish_vals, scenario_3_coral, 100, colour_scheme(3, :), 'filled')
s = scatter(initial_starfish_vals, scenario_4_coral, 'x');
s.MarkerFaceColor = colour_scheme(4, :);
s.MarkerEdgeColor = colour_scheme(4, :);
s.SizeData = 100;
s.LineWidth = 2;
set(gca, 'FontSize', ticks_FS);
xlabel('Initial no. of adult starfish', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total coral cover', 'Interpreter', 'Latex', ...
    'Fontsize', axis_FS)
title('\qquad\qquad\qquad\qquad\qquad Total coral cover on GBR after 100 years for varying initial no. of adult starfish', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
legend('No control', '100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
    '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs',  ...
    'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)