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
t_end = 200;                     % time in years

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

% Connectivity matrices from Bode et al. (2012)
params.omega_c = psurv_d02_1122_P7;             % coral
params.omega_s = psurv_d02_1122_P7;             % starfish

% Latitude and longitude for starfish larval calculation
params.lon = lat;
params.lat = lon;

% Using metapopulation model equation for larval dispersal
dispersal_eq = 1;


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


%% SOLVE ==================================================================
% SETUP -------------------------------------------------------------------
% Setup fecundity parameter arrays
r_c_vals = 0:0.1:5;            % coral larvae reproduction rate
r_s_vals = 0:100:5000;         % starfish larvae reproduction rate

% Initialise multi-dimensional arrays to store solution
coral = zeros(length(r_c_vals), length(r_s_vals), t_end+1);
starfish_age2 = zeros(length(r_c_vals), length(r_s_vals), t_end+1);
starfish_age1 = zeros(length(r_c_vals), length(r_s_vals), t_end+1);
starfish_age0 = zeros(length(r_c_vals), length(r_s_vals), t_end+1);


% LOOP OVER PARAMETERS ----------------------------------------------------
for i = 1:length(r_c_vals)
    for j = 1:length(r_s_vals)
        % CONTROL SCENARIO 0 ----------------------------------------------
        % SOLVE 
        % Choose parameter values
        params.r_c = r_c_vals(i);
        params.r_s = r_s_vals(j);
        
        % Solve using function which runs simulations
        [t_vec_s0, C_y_f_s0, N_y_2_s0, N_y_1_s0, N_y_0_s0, ~] = ...
            simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s0, dispersal_eq);

        % CALCULATIONS 
        % Calculate coral cover and cots over time
        coral_s0 = sum(C_y_f_s0, 1);
        starfish_age2_s0 = sum(N_y_2_s0, 1);
        starfish_age1_s0 = sum(N_y_1_s0, 1);
        starfish_age0_s0 = sum(N_y_0_s0, 1);

        % Calculate coral cover in initiation box over time
        [coral_box_s0, starfish_age2_box_s0, starfish_age1_box_s0, starfish_age0_box_s0] = ...
            calculate_population_box(t_end, C_y_f_s0, N_y_2_s0, N_y_1_s0, N_y_0_s0, num_reefs, lon, lat);
        
        % SAVE
        coral(i, j, :) = coral_s0;
        starfish_age2(i, j, :) = starfish_age2_s0;
        starfish_age1(i, j, :) = starfish_age1_s0;
        starfish_age0(i, j, :) = starfish_age0_s0;
    end
end



%% PLOTS ==================================================================

% Fontsizes for plotting
axis_FS = 14;
title_FS = 15;
legend_FS = 12;
ticks_FS = 12;
colorbar_FS = 14;

for i = 1:t_end+1
    % Plot coral at final timestep
    figure(10), clf, hold on
    sp1 = subplot(2, 1, 1); hold on
    p1 = surf(r_s_vals, r_c_vals, coral(:, :, i));
    set(p1, 'EdgeColor', 'none');
    shading interp
    view(2)
    c = colorbar;
    colormap(sp1, viridis)
    set(gca, 'FontSize', ticks_FS);
    xlabel('$\lambda_s$', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
    ylabel('$\lambda_c$', 'Interpreter', 'Latex', ...
        'Fontsize', axis_FS)
    title(['Total coral cover on GBR after ', num2str(i-1), ' years'], ...
        'Interpreter', 'Latex', 'Fontsize', title_FS)
    c.Label.String = 'Total coral cover';
    c.Label.Interpreter = 'Latex';
    c.Label.FontSize = colorbar_FS;
%     caxis([0 2750])

    % Plot starfish at final timestep
    sp2 = subplot(2, 1, 2); hold on
    p2 = surf(r_s_vals(16:31), r_c_vals, starfish_age2(:, 16:31, i));
%     p2 = surf(r_s_vals, r_c_vals, starfish_age2(:, :, i));
    set(p2, 'EdgeColor', 'none');
    shading interp
    view(2)
    c = colorbar;
    colormap(sp2, magma)
    set(gca, 'FontSize', ticks_FS);
    xlabel('$\lambda_s$', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
    ylabel('$\lambda_c$', 'Interpreter', 'Latex', ...
        'Fontsize', axis_FS)
    title(['Total age 2+ starfish on GBR after ', num2str(i-1), ' years'], ...
        'Interpreter', 'Latex', 'Fontsize', title_FS)
    c.Label.String = 'No. of age 2+ starfish';
    c.Label.Interpreter = 'Latex';
    c.Label.FontSize = colorbar_FS;
%     caxis([0 12e5])
end 

saveas(gcf, 'Plots/05_sensitivity_analysis/lambda_parameters.png')
