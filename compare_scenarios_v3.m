%% SETUP
clear all

% MODEL ===================================================================
% DATA --------------------------------------------------------------------
% Load in the connectivity matrices
load IdentifyKeySources/ConnectivityMatrices_Model_A_2002_P7

% Load in Australian map outline
load DataSources/AustOutline

% Load in latitude longitude data
load IdentifyKeySources/original_centroids

% Match the reefs to the areas
[num_reefs, lon, lat, reef_area, omega] = match_reefs_to_areas(lg, lt, ...
    psurv_d02_1122_P7);

% Convert reef area from hectares to sqkm
reef_area = reef_area / 100;

% Clear the rest of the connectivity matrices that we don't need
clear psurv*


% PARAMETERS --------------------------------------------------------------
% How long do we want to run the simulation for
t_end = 100;                    % time in years

% Parameters struct: store all parameters in a struct
% Reef areas taken from Ryu's data
params.K_f = reef_area;         % carrying capacity of fast-growing coral

% Estimated by Morello et al. (2014)
params.p_tilde = 0.258;         % effect of fast-growing coral on COTS     
params.M_cots = 2.56;           % natural mortality of COTS 
params.p_1_f = ...
    (0.129/2500)*params.K_f;    % effect of COTS on fast-growing coral

% Known or arbitrarily chosen by Morello et al. (2014)
params.r_f = 0.5;               % intrinsic growth rate of fast-growing coral
params.p_2_f = ...
    (10/2500)*params.K_f;       % effect of COTS on fast-growing coral

% Arbitrarily chosen by me
params.r_c = 0.1;               % coral larvae reproduction rate
params.r_s = 5000;              % starfish larvae reproduction rate

% Call function to calculate percentage of age 1 COTS that reproduce, 
% calculated from Lucas (1984) and Babcock et al. (2016)
params.mu_s = calculate_cots_age1_reproduction();     

% Connectivity matrices from Bode et al. (2012)
V_f = 0.9;                      % coral larval survival rate
V_s = 0.25;                      % starfish larval survival rate
params.omega_c = V_f*omega;     % coral larval dispersal
params.omega_s = V_s*omega;     % starfish larval dispersal

% Latitude and longitude for starfish larval calculation
params.lon = lon;
params.lat = lat;

% Using metapopulation model equation for larval dispersal
dispersal_eq = 1;


% INITIAL SYSTEM STATE ----------------------------------------------------
% CORAL
% Percentage of fast-growing coral = 50% of reef area everywhere
initial_coral = 0.5;
initial_state.C_0_f = initial_coral * params.K_f;

% Calculate carrying capacity of coral on entire GBR
coral_GBR_cc = sum(reef_area);

% Calculate carrying capacity of coral at initiation box
coral_box_cc = 0;
for i = 1:num_reefs
    if (lat(i) > -17 && lat(i) < -14.75) && (lon(i) > 145 && lon(i) < 147)
        coral_box_cc = coral_box_cc + reef_area(i);
    end
end

% STARFISH
% Initialise array for number of age 2+ COTS
initial_state.N_0_2 = zeros(num_reefs, 1);

% Look for reefs within the initiation box, and put some starfish there
for i = 1:num_reefs
    if (lat(i) > -17 && lat(i) < -14.75) && (lon(i) > 145 && lon(i) < 147)
        initial_state.N_0_2(i) = 100;
    end
end

% Initialise age 0 and age 1 COTS based on modified Morello initial conditions
initial_state.N_0_1 = initial_state.N_0_2 .* exp(params.M_cots * (1 ...
    - params.p_tilde * (initial_state.C_0_f ./ (1 + initial_state.C_0_f))));
initial_state.N_0_0 = initial_state.N_0_2 .* exp(2 * params.M_cots * (1 ...
    - params.p_tilde * (initial_state.C_0_f ./ (1 + initial_state.C_0_f))));



%% SCENARIO 0: No control, but simulation run with same conditions

% CONTROL EFFORT ----------------------------------------------------------
% No control
control_effort_s0 = 0;

% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s0, C_y_f_s0, N_y_2_s0, N_y_1_s0, N_y_0_s0, tau_ratio_s0] = ...
    simulate_reefs_v3(num_reefs, t_end, params, initial_state, control_effort_s0, dispersal_eq);

% CALCULATIONS ------------------------------------------------------------
% Calculate coral cover and cots over time
coral_s0 = sum(C_y_f_s0, 1);
starfish_age2_s0 = sum(N_y_2_s0, 1);
starfish_age1_s0 = sum(N_y_1_s0, 1);
starfish_age0_s0 = sum(N_y_0_s0, 1);

% Calculate coral cover in initiation box over time
[coral_box_s0, starfish_age2_box_s0, starfish_age1_box_s0, starfish_age0_box_s0] = ...
    calculate_population_box(t_end, C_y_f_s0, N_y_2_s0, N_y_1_s0, N_y_0_s0, num_reefs, lon, lat);



%% SCENARIO 1: Control at initiation box

% CONTROL EFFORT ----------------------------------------------------------
% Cull all the starfish only in the initiation box for all 100 years
control_effort_s1 = zeros(num_reefs, t_end);
for i = 1:num_reefs
    if (lat(i) > -17 && lat(i) < -14.75) && (lon(i) > 145 && lon(i) < 147)
        control_effort_s1(i, :) = 1;
    end
end

% Calculate number of reefs controlled, area controlled and budget
num_reefs_control_s1 = nnz(control_effort_s1(:, 1))
area_control_s1 = sum(reef_area .* control_effort_s1(:, 1) * 1, 'all')
budget_s1 = sum(repmat(reef_area, 1, 100) .* control_effort_s1, 'all')


% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s1, C_y_f_s1, N_y_2_s1, N_y_1_s1, N_y_0_s1, tau_ratio_s1] = ...
    simulate_reefs_v3(num_reefs, t_end, params, initial_state, control_effort_s1, dispersal_eq);


% CALCULATIONS ------------------------------------------------------------
% Calculate coral cover and cots over time
coral_s1 = sum(C_y_f_s1, 1);
starfish_age2_s1 = sum(N_y_2_s1, 1);
starfish_age1_s1 = sum(N_y_1_s1, 1);
starfish_age0_s1 = sum(N_y_0_s1, 1);

% Calculate coral cover in initiation box over time
[coral_box_s1, starfish_age2_box_s1, starfish_age1_box_s1, starfish_age0_box_s1] = ...
    calculate_population_box(t_end, C_y_f_s1, N_y_2_s1, N_y_1_s1, N_y_0_s1, num_reefs, lon, lat);



%% SCENARIO 2: Control at twice the area but half the effort 

% CONTROL EFFORT ----------------------------------------------------------
% Cull starfish in twice the area at half the effort
control_effort_s2 = zeros(num_reefs, t_end);
for i = 1:num_reefs
    if (lat(i) > -18.045 && lat(i) < -13.98) && (lon(i) > 142 && lon(i) < 147.5)
        control_effort_s2(i, :) = 0.5;
    end
end

% Calculate number of reefs controlled, area controlled and budget
num_reefs_control_s2 = nnz(control_effort_s2(:, 1))
area_control_s2 = sum(reef_area .* control_effort_s2(:, 1) * 2, 'all')
budget_s2 = sum(repmat(reef_area, 1, 100) .* control_effort_s2, 'all') 


% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s2, C_y_f_s2, N_y_2_s2, N_y_1_s2, N_y_0_s2, tau_ratio_s2] = ...
    simulate_reefs_v3(num_reefs, t_end, params, initial_state, control_effort_s2, dispersal_eq);


% CALCULATIONS ------------------------------------------------------------
% Calculate coral cover and cots over time
coral_s2 = sum(C_y_f_s2, 1);
starfish_age2_s2 = sum(N_y_2_s2, 1);
starfish_age1_s2 = sum(N_y_1_s2, 1);
starfish_age0_s2 = sum(N_y_0_s2, 1);

% Calculate coral cover in initiation box over time
[coral_box_s2, starfish_age2_box_s2, starfish_age1_box_s2, starfish_age0_box_s2] = ...
    calculate_population_box(t_end, C_y_f_s2, N_y_2_s2, N_y_1_s2, N_y_0_s2, num_reefs, lon, lat);



%% SCENARIO 3: Control quadruple the area but same total effort

% CONTROL EFFORT ----------------------------------------------------------
% Cull starfish in four times the area and a quarter of effort at each reef
control_effort_s3 = zeros(num_reefs, t_end);
for i = 1:num_reefs
    if (lat(i) > -19.73 && lat(i) < -13) && (lon(i) > 142 && lon(i) < 152)
        control_effort_s3(i, :) = 0.25;
    end
end

% Calculate number of reefs controlled, area controlled and budget
num_reefs_control_s3 = nnz(control_effort_s3(:, 1))
area_control_s3 = sum(reef_area .* control_effort_s3(:, 1) * 4, 'all')
budget_s3 = sum(repmat(reef_area, 1, 100) .* control_effort_s3, 'all')


% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s3, C_y_f_s3, N_y_2_s3, N_y_1_s3, N_y_0_s3, tau_ratio_s3] = ...
    simulate_reefs_v3(num_reefs, t_end, params, initial_state, control_effort_s3, dispersal_eq);


% CALCULATIONS ------------------------------------------------------------
% Calculate coral cover and cots over time
coral_s3 = sum(C_y_f_s3, 1);
starfish_age2_s3 = sum(N_y_2_s3, 1);
starfish_age1_s3 = sum(N_y_1_s3, 1);
starfish_age0_s3 = sum(N_y_0_s3, 1);

% Calculate coral cover in initiation box over time
[coral_box_s3, starfish_age2_box_s3, starfish_age1_box_s3, starfish_age0_box_s3] = ...
    calculate_population_box(t_end, C_y_f_s3, N_y_2_s3, N_y_1_s3, N_y_0_s3, num_reefs, lon, lat);



%% SCENARIO 4: Control entire reef with the same budget as others

% CONTROL EFFORT ----------------------------------------------------------
% Cull starfish everywhere based on total budget available
effort_per_year = budget_s1 / t_end;
effort_per_reef = effort_per_year / coral_GBR_cc
control_effort_s4 = effort_per_reef * ones(num_reefs, t_end);

% Calculate number of reefs controlled, area controlled and budget
num_reefs_control_s4 = nnz(control_effort_s4(:, 1))
area_control_s4 = sum(reef_area)
budget_s4 = sum(repmat(reef_area, 1, 100) .* control_effort_s4, 'all')


% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s4, C_y_f_s4, N_y_2_s4, N_y_1_s4, N_y_0_s4, tau_ratio_s4] = ...
    simulate_reefs_v3(num_reefs, t_end, params, initial_state, control_effort_s4, dispersal_eq);


% CALCULATIONS ------------------------------------------------------------
% Calculate coral cover and cots over time
coral_s4 = sum(C_y_f_s4, 1);
starfish_age2_s4 = sum(N_y_2_s4, 1);
starfish_age1_s4 = sum(N_y_1_s4, 1);
starfish_age0_s4 = sum(N_y_0_s4, 1);

% Calculate coral cover in initiation box over time
[coral_box_s4, starfish_age2_box_s4, starfish_age1_box_s4, starfish_age0_box_s4] = ...
    calculate_population_box(t_end, C_y_f_s4, N_y_2_s4, N_y_1_s4, N_y_0_s4, num_reefs, lon, lat);



%% OTHER SCENARIO: Control a ridiculous amount to see coral cover not die

% CONTROL EFFORT ----------------------------------------------------------
% Cull all the adult starfish everywhere
control_effort_other = ones(num_reefs, t_end);

% Calculate area controlled and budget 
area_control_other = sum(reef_area)
budget_other = sum(repmat(reef_area, 1, 100) .* control_effort_other, 'all')

% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_other, C_y_f_other, N_y_2_other, N_y_1_other, N_y_0_other, tau_ratio_other] = ...
    simulate_reefs_v3(num_reefs, t_end, params, initial_state, control_effort_other, dispersal_eq);

% CALCULATIONS ------------------------------------------------------------
% Calculate coral cover and cots over time
coral_other = sum(C_y_f_other, 1);
starfish_age2_other = sum(N_y_2_other, 1);
starfish_age1_other = sum(N_y_1_other, 1);
starfish_age0_other = sum(N_y_0_other, 1);

% Calculate coral cover in initiation box over time
[coral_box_other, starfish_age2_box_other, starfish_age1_box_other, starfish_age0_box_other] = ...
    calculate_population_box(t_end, C_y_f_other, N_y_2_other, N_y_1_other, N_y_0_other, num_reefs, lon, lat);



%% COMPARE: Compare total coral cover at end of each scenario

% Bar graph calculations --------------------------------------------------
% Get the coral cover at the end of 100 years
coral_end_s0 = C_y_f_s0(:, end) ./ reef_area;
coral_end_s1 = C_y_f_s1(:, end) ./ reef_area;
coral_end_s2 = C_y_f_s2(:, end) ./ reef_area;
coral_end_s3 = C_y_f_s3(:, end) ./ reef_area;
coral_end_s4 = C_y_f_s4(:, end) ./ reef_area;

% Combine into matrix
coral_end = [coral_end_s0 coral_end_s1 coral_end_s2 coral_end_s3 coral_end_s4];

% Normalise
coral_end_norm = normalize(coral_end, 'center');

% Initialise matrix to store comparisons
coral_compare_all = zeros(size(coral_end, 2), 5);

% Loop over each column in matrix and count number of reefs
for i = 1:size(coral_end, 2)
    % Count the number of reefs with less than 10% coral 
    coral_compare_all(i, 1) = sum(coral_end(:, i) <= 0.1);
    
    % Count the number of reefs with between 10% and 30% coral
    coral_compare_all(i, 2) = sum(coral_end(:, i) > 0.1 & coral_end(:, i) <= 0.3);
    
    % Count the number of reefs with between 30% and 50% coral
    coral_compare_all(i, 3) = sum(coral_end(:, i) > 0.3 & coral_end(:, i) <= 0.5);
    
    % Count the number of reefs with between 50% and 75% coral
    coral_compare_all(i, 4) = sum(coral_end(:, i) > 0.5 & coral_end(:, i) <= 0.75);
    
    % Count the number of reefs with more than 75% coral
    coral_compare_all(i, 5) = sum(coral_end(:, i) > 0.75);
end

% Count how many extra reefs in each category compared to no coral
coral_compare = coral_compare_all(2:5, :) - coral_compare_all(1, :);


% Histogram calculations --------------------------------------------------
% Combine coral at end of scenarios into matrix
coral_end_control = [coral_end_s1 coral_end_s2 coral_end_s3 coral_end_s4];

% Initialise matrix
coral_end_diff = zeros(num_reefs, size(coral_end_control, 2));

% Find difference in coral cover
for i = 1:size(coral_end_control, 2)
    coral_end_diff(:, i) = coral_end_control(:, i) - coral_end_s0;
end

% Normalise 
coral_end_diff_norm = normalize(coral_end_diff, 'center');



%% COMPARE: Compare the control effort in each scenario

% Get the control effort values for only one timestep
control_1_plot = control_effort_s1(:, 1);
control_2_plot = control_effort_s2(:, 1);
control_3_plot = control_effort_s3(:, 1);
control_4_plot = control_effort_s4(:, 1);

% Combine into matrix
control_plot = [control_1_plot control_2_plot control_3_plot control_4_plot];



%% CALCULATE: Histogram calculations for latitudinal spread



%% PLOT FORMATTING ========================================================
% Fontsizes for plotting
axis_FS = 15;
title_FS = 16;
legend_FS = 13;
ticks_FS = 12;

% Colours for plotting
colour_scheme = cbrewer('div', 'RdBu', 4);
colour_scheme_brown = cbrewer('div', 'BrBG', 9);
colour_scheme_purple = cbrewer('div', 'PRGn', 9);

% Define viridis and magma palette colours 
viridis_palette_4 = [253, 231, 37; 53, 183, 121; 49, 104, 142; 68, 1, 84];
viridis_palette_4 = viridis_palette_4/255;
viridis_palette_5 = [253, 231, 37; 94, 201, 98; 33, 145, 140; 59, 82, 139; 68, 1, 84];
viridis_palette_5 = viridis_palette_5/255;
magma_palette = [252, 253, 191; 252, 137, 97; 183, 55, 121; 81, 18, 124; 0, 0, 4];
magma_palette = magma_palette/255;



%% CORAL ==================================================================
% GBR total coral cover ---------------------------------------------------
figure(1), clf, hold on
yline(coral_GBR_cc, '--', 'Linewidth', 2, 'Color', [0.5 0.5 0.5])
plot(t_vec_s0, coral_s0, 'Linewidth', 2, 'Color', 'black')
plot(t_vec_other, coral_other, 'Linewidth', 2, 'Color', colour_scheme_brown(1, :))
plot(t_vec_s1, coral_s1, 'Linewidth', 2, 'Color', colour_scheme(1, :))
plot(t_vec_s2, coral_s2, 'Linewidth', 2, 'Color', colour_scheme(2, :))
plot(t_vec_s3, coral_s3, 'Linewidth', 2, 'Color', colour_scheme(3, :))
plot(t_vec_s4, coral_s4, 'Linewidth', 2, 'Color', colour_scheme(4, :))
set(gca, 'FontSize', ticks_FS);
ylim([0 16000])
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total coral cover ($km^2$)', 'Interpreter', 'Latex', ...
    'Fontsize', axis_FS)
title('Total coral cover on GBR', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
legend('Coral carrying capacity', 'No control', '100\% effort over 14240 $km^2$', ...
    '100\% effort over 1737.1 $km^2$', '50\% effort over 3474.2 $km^2$', ...
    '25\% effort over 6948.4 $km^2$', '12.2\% effort over 14240 $km^2$',  ...
    'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% GBR total coral cover (additional) --------------------------------------
figure(2), clf, hold on
plot(t_vec_s1, coral_s1-coral_s0, 'Linewidth', 2, 'Color', colour_scheme(1, :))
plot(t_vec_s2, coral_s2-coral_s0, 'Linewidth', 2, 'Color', colour_scheme(2, :))
plot(t_vec_s3, coral_s3-coral_s0, 'Linewidth', 2, 'Color', colour_scheme(3, :))
plot(t_vec_s4, coral_s4-coral_s0, 'Linewidth', 2, 'Color', colour_scheme(4, :))
set(gca, 'FontSize', ticks_FS);
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total coral cover increase ($km^2$)', 'Interpreter', 'Latex', ...
    'Fontsize', axis_FS)
title('Total coral cover increase on GBR compared to no control', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
legend('100\% effort over 1737.1 $km^2$', '50\% effort over 3474.2 $km^2$', ...
    '25\% effort over 6948.4 $km^2$', '12.2\% effort over 14240 $km^2$',  ...
    'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% Initiation box total coral cover ----------------------------------------
figure(3), clf, hold on
yline(coral_box_cc, '--', 'Linewidth', 2, 'Color', [0.5 0.5 0.5])
plot(t_vec_s0, coral_box_s0, 'Linewidth', 2, 'Color', 'black')
plot(t_vec_s1, coral_box_s1, 'Linewidth', 2, 'Color', colour_scheme(1, :))
plot(t_vec_s2, coral_box_s2, 'Linewidth', 2, 'Color', colour_scheme(2, :))
plot(t_vec_s3, coral_box_s3, 'Linewidth', 2, 'Color', colour_scheme(3, :))
plot(t_vec_s4, coral_box_s4, 'Linewidth', 2, 'Color', colour_scheme(4, :))
set(gca, 'FontSize', ticks_FS);
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total coral cover ($km^2$)', 'Interpreter', 'Latex', ...
    'Fontsize', axis_FS)
title('Total coral cover at initiation box', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
legend('Coral carrying capacity', 'No control', '100\% effort over 1737.1 $km^2$', ...
    '50\% effort over 3474.2 $km^2$', '25\% effort over 6948.4 $km^2$', ...
    '12.2\% effort over 14240 $km^2$', 'Location', 'NorthEastOutside', ...
    'Interpreter', 'Latex', 'Fontsize', legend_FS)

% Initiation box total coral cover (additional) ---------------------------
figure(4), clf, hold on
plot(t_vec_s1, coral_box_s1-coral_box_s0, 'Linewidth', 2, 'Color', colour_scheme(1, :))
plot(t_vec_s2, coral_box_s2-coral_box_s0, 'Linewidth', 2, 'Color', colour_scheme(2, :))
plot(t_vec_s3, coral_box_s3-coral_box_s0, 'Linewidth', 2, 'Color', colour_scheme(3, :))
plot(t_vec_s4, coral_box_s4-coral_box_s0, 'Linewidth', 2, 'Color', colour_scheme(4, :))
set(gca, 'FontSize', ticks_FS);
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total coral cover increase ($km^2$)', 'Interpreter', 'Latex', ...
    'Fontsize', axis_FS)
title('\qquad\qquad Total coral cover increase at initiation box compared to no control', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
legend('100\% effort over 1737.1 $km^2$', '50\% effort over 3474.2 $km^2$', ...
    '25\% effort over 6948.4 $km^2$', '12.2\% effort over 14240 $km^2$',  ...
    'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)



%% STARFISH ===============================================================
% GBR total age 2+ starfish -----------------------------------------------
figure(5), clf, hold on
plot(t_vec_s0, starfish_age2_s0, 'Linewidth', 2, 'Color', 'Black')
plot(t_vec_s1, starfish_age2_s1, 'Linewidth', 2, 'Color', colour_scheme(1, :))
plot(t_vec_s2, starfish_age2_s2, 'Linewidth', 2, 'Color', colour_scheme(2, :))
plot(t_vec_s3, starfish_age2_s3, 'Linewidth', 2, 'Color', colour_scheme(3, :))
plot(t_vec_s4, starfish_age2_s4, 'Linewidth', 2, 'Color', colour_scheme(4, :))
set(gca, 'FontSize', ticks_FS)
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('No. of adult starfish', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title('Total adult starfish population on GBR', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
legend('No control', '100\% effort over 1737.1 $km^2$', ...
    '50\% effort over 3474.2 $km^2$', '25\% effort over 6948.4 $km^2$', ...
    '12.2\% effort over 14240 $km^2$', 'Location', 'NorthEastOutside', ...
    'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % GBR total age 2+ starfish (log scale) -----------------------------------
% figure(5), clf, hold on
% plot(t_vec_s0, starfish_age2_s0, 'Linewidth', 2, 'Color', 'Black')
% plot(t_vec_s1, starfish_age2_s1, 'Linewidth', 2, 'Color', colour_scheme(1, :))
% plot(t_vec_s2, starfish_age2_s2, 'Linewidth', 2, 'Color', colour_scheme(2, :))
% plot(t_vec_s3, starfish_age2_s3, 'Linewidth', 2, 'Color', colour_scheme(3, :))
% plot(t_vec_s4, starfish_age2_s4, 'Linewidth', 2, 'Color', colour_scheme(4, :))
% set(gca, 'YScale', 'log')
% set(gca, 'FontSize', ticks_FS)
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('Total no. of adult starfish (log scale)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title('Total adult starfish population on GBR', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('No control', '100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% GBR decrease in age 2+ starfish -----------------------------------------
figure(6), clf, hold on
plot(t_vec_s1, starfish_age2_s0-starfish_age2_s1, 'Linewidth', 2, 'Color', colour_scheme(1, :))
plot(t_vec_s2, starfish_age2_s0-starfish_age2_s2, 'Linewidth', 2, 'Color', colour_scheme(2, :))
plot(t_vec_s3, starfish_age2_s0-starfish_age2_s3, 'Linewidth', 2, 'Color', colour_scheme(3, :))
plot(t_vec_s4, starfish_age2_s0-starfish_age2_s4, 'Linewidth', 2, 'Color', colour_scheme(4, :))
% set(gca, 'YScale', 'log')
set(gca, 'FontSize', ticks_FS)
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total adult starfish decrease', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title('\qquad\qquad\qquad\qquad\qquad\qquad Total decrease in adult starfish population on GBR compared to no control', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
legend('100\% effort over 1737.1 $km^2$', '50\% effort over 3474.2 $km^2$', ...
    '25\% effort over 6948.4 $km^2$', '12.2\% effort over 14240 $km^2$',  ...
    'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % GBR total age 2+ starfish (control scenarios only) ----------------------
% figure(6), clf, hold on, grid on
% plot(t_vec_s1, starfish_age2_s1, 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980])
% plot(t_vec_s2, starfish_age2_s2, 'Linewidth', 2, 'Color', [0.9290 0.6940 0.1250])
% plot(t_vec_s3, starfish_age2_s3, 'Linewidth', 2, 'Color', [0.4940 0.1840 0.5560])
% plot(t_vec_s4, starfish_age2_s4, 'Linewidth', 2, 'Color', [0.4660 0.6740 0.1880])
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 2+ starfish', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title('Total adult starfish population on GBR', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% Initiation box total age 2+ starfish ------------------------------------
figure(7), clf, hold on
plot(t_vec_s0, starfish_age2_box_s0, 'Linewidth', 2, 'Color', 'black')
plot(t_vec_s1, starfish_age2_box_s1, 'Linewidth', 2, 'Color', colour_scheme(1, :))
plot(t_vec_s2, starfish_age2_box_s2, 'Linewidth', 2, 'Color', colour_scheme(2, :))
plot(t_vec_s3, starfish_age2_box_s3, 'Linewidth', 2, 'Color', colour_scheme(3, :))
plot(t_vec_s4, starfish_age2_box_s4, 'Linewidth', 2, 'Color', colour_scheme(4, :))
set(gca, 'FontSize', ticks_FS)
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('No. of adult starfish', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title('Total adult starfish population at initiation box', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
legend('No control', '100\% effort over 1737.1 $km^2$', ...
    '50\% effort over 3474.2 $km^2$', '25\% effort over 6948.4 $km^2$', ...
    '12.2\% effort over 14240 $km^2$', 'Location', 'NorthEastOutside', ...
    'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Initiation box total age 2+ starfish (log scale) ------------------------
% figure(7), clf, hold on
% plot(t_vec_s0, starfish_age2_box_s0, 'Linewidth', 2, 'Color', 'Black')
% plot(t_vec_s1, starfish_age2_box_s1, 'Linewidth', 2, 'Color', colour_scheme(1, :))
% plot(t_vec_s2, starfish_age2_box_s2, 'Linewidth', 2, 'Color', colour_scheme(2, :))
% plot(t_vec_s3, starfish_age2_box_s3, 'Linewidth', 2, 'Color', colour_scheme(3, :))
% plot(t_vec_s4, starfish_age2_box_s4, 'Linewidth', 2, 'Color', colour_scheme(4, :))
% set(gca, 'FontSize', ticks_FS)
% % set(gca, 'YGrid', 'on')
% set(gca, 'YScale', 'log')
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('Total no. of adult starfish', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title('Total adult starfish population at initiation box', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('No control', '100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% Initiation box decrease in age 2+ starfish ------------------------------
figure(8), clf, hold on
plot(t_vec_s1, starfish_age2_box_s0-starfish_age2_box_s1, 'Linewidth', 2, 'Color', colour_scheme(1, :))
plot(t_vec_s2, starfish_age2_box_s0-starfish_age2_box_s2, 'Linewidth', 2, 'Color', colour_scheme(2, :))
plot(t_vec_s3, starfish_age2_box_s0-starfish_age2_box_s3, 'Linewidth', 2, 'Color', colour_scheme(3, :))
plot(t_vec_s4, starfish_age2_box_s0-starfish_age2_box_s4, 'Linewidth', 2, 'Color', colour_scheme(4, :))
% set(gca, 'YScale', 'log')
set(gca, 'FontSize', ticks_FS)
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total adult starfish decrease', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title('\quad\qquad\qquad\qquad\qquad\qquad Total decrease in adult starfish population at initiation box compared to no control', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
legend('100\% effort over 1737.1 $km^2$', '50\% effort over 3474.2 $km^2$', ...
    '25\% effort over 6948.4 $km^2$', '12.2\% effort over 14240 $km^2$',  ...
    'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Initiation box total age 2+ starfish (additional) -----------------------
% figure(9), clf, hold on, grid on
% plot(t_vec_s0, starfish_age2_box_s0-starfish_age2_box_s0, 'Linewidth', 2)
% plot(t_vec_s1, starfish_age2_box_s1-starfish_age2_box_s0, 'Linewidth', 2)
% plot(t_vec_s2, starfish_age2_box_s2-starfish_age2_box_s0, 'Linewidth', 2)
% plot(t_vec_s3, starfish_age2_box_s3-starfish_age2_box_s0, '--', 'Linewidth', 2)
% plot(t_vec_s4, starfish_age2_box_s4-starfish_age2_box_s0, 'Linewidth', 2)
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 2+ starfish', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title('\qquad\qquad\qquad Additional total adult starfish population in initiation box', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('No control', '100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Initiation box total age 2+ starfish (control scenarios only) -----------
% figure(10), clf, hold on, grid on
% plot(t_vec_s1, starfish_age2_box_s1, 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980])
% plot(t_vec_s2, starfish_age2_box_s2, 'Linewidth', 2, 'Color', [0.9290 0.6940 0.1250])
% plot(t_vec_s3, starfish_age2_box_s3, 'Linewidth', 2, 'Color', [0.4940 0.1840 0.5560])
% plot(t_vec_s4, starfish_age2_box_s4, 'Linewidth', 2, 'Color', [0.4660 0.6740 0.1880])
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 2+ starfish', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title('Total adult starfish population in initiation box', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Initiation box total age 1 starfish (control scenarios only) ------------
% figure(11), clf, hold on, grid on
% plot(t_vec_s1, starfish_age1_box_s1, 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980])
% plot(t_vec_s2, starfish_age1_box_s2, 'Linewidth', 2, 'Color', [0.9290 0.6940 0.1250])
% plot(t_vec_s3, starfish_age1_box_s3, 'Linewidth', 2, 'Color', [0.4940 0.1840 0.5560])
% plot(t_vec_s4, starfish_age1_box_s4, 'Linewidth', 2, 'Color', [0.4660 0.6740 0.1880])
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 1 starfish', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title('Total age 1 starfish population in initiation box', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Initiation box total age 0 starfish (control scenarios only) ------------
% figure(12), clf, hold on, grid on
% plot(t_vec_s1, starfish_age0_box_s1, 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980])
% plot(t_vec_s2, starfish_age0_box_s2, 'Linewidth', 2, 'Color', [0.9290 0.6940 0.1250])
% plot(t_vec_s3, starfish_age0_box_s3, 'Linewidth', 2, 'Color', [0.4940 0.1840 0.5560])
% plot(t_vec_s4, starfish_age0_box_s4, 'Linewidth', 2, 'Color', [0.4660 0.6740 0.1880])
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 0 starfish', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title('Total age 0 starfish population in initiation box', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Initiation box total starfish (control scenarios only) ------------------
% figure(13), clf, hold on, grid on
% plot(t_vec_s1, starfish_age0_box_s1+starfish_age1_box_s1+starfish_age2_box_s1, ...
%     'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980])
% plot(t_vec_s2, starfish_age0_box_s2+starfish_age1_box_s2+starfish_age2_box_s2, ...
%     'Linewidth', 2, 'Color', [0.9290 0.6940 0.1250])
% plot(t_vec_s3, starfish_age0_box_s3+starfish_age1_box_s3+starfish_age2_box_s3, ...
%     'Linewidth', 2, 'Color', [0.4940 0.1840 0.5560])
% plot(t_vec_s4, starfish_age0_box_s4+starfish_age1_box_s4+starfish_age2_box_s4, ...
%     'Linewidth', 2, 'Color', [0.4660 0.6740 0.1880])
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of starfish', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title('Total starfish population in initiation box', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Starfish larval dispersal ratio -----------------------------------------
% figure(14), clf, hold on, grid on
% plot(t_vec_s0, tau_ratio_s0, 'Linewidth', 2)
% plot(t_vec_s1, tau_ratio_s1, 'Linewidth', 2)
% plot(t_vec_s2, tau_ratio_s2, 'Linewidth', 2)
% plot(t_vec_s3, tau_ratio_s3, 'Linewidth', 2)
% plot(t_vec_s4, tau_ratio_s4, 'Linewidth', 2)
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('\% of starfish larvae', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title('\qquad\qquad\qquad Percentage of age 0 starfish arriving at initiation box, born at initiation box', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('No control', '100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 rgueefs', ...
%     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)



%% HEATMAPS ===============================================================
% % Coral heatmap -----------------------------------------------------------
% f1 = figure(15); clf
% imagesc(C_y_f_s0)
% title('Coral cover on GBR', 'Interpreter', 'Latex', 'Fontsize', title_FS)
% colorbar
% colormap(f1, viridis)

% % Starfish heatmap --------------------------------------------------------
% f2 = figure(16); clf
% imagesc(N_y_2_s0)
% title('Adult starfish on GBR', 'Interpreter', 'Latex', 'Fontsize', title_FS)
% colorbar
% colormap(f2, magma)



%% NUMBER OF REEFS WITH STARFISH ==========================================
% % Calculate number of reefs for each age class
% for i = 1:t_end+1
%     cots_age2_reefs(i) = nnz(N_y_2_s0(:, i));
%     cots_age1_reefs(i) = nnz(N_y_1_s0(:, i));
%     cots_age0_reefs(i) = nnz(N_y_0_s0(:, i));
% end
%
% % Number of reefs with starfish -------------------------------------------
% figure(17), clf, hold on
% plot(t_vec_s0, cots_age2_reefs, 'Linewidth', 2)
% plot(t_vec_s0, cots_age1_reefs, 'Linewidth', 2)
% plot(t_vec_s0, cots_age0_reefs, 'Linewidth', 2)
% legend('Age 2+', 'Age 1', 'Age 0')
% xlabel('Time')
% ylabel('Number of reefs')
% title('Number of reefs with starfish over time')
% xlim([0 20])



%% SINGLE REEF ============================================================
% % Pick reef to plot 
% reef = 184;
% 
% % Coral cover at one reef -------------------------------------------------
% figure(18), clf, hold on, grid on
% plot(t_vec_s0, C_y_f_s0(reef, :), 'Linewidth', 2, 'Color', 'black')
% plot(t_vec_s1, C_y_f_s1(reef, :), 'Linewidth', 2, 'Color', colour_scheme(1, :))
% plot(t_vec_s2, C_y_f_s2(reef, :), 'Linewidth', 2, 'Color', colour_scheme(2, :))
% plot(t_vec_s3, C_y_f_s3(reef, :), 'Linewidth', 2, 'Color', colour_scheme(3, :))
% plot(t_vec_s4, C_y_f_s4(reef, :), 'Linewidth', 2, 'Color', colour_scheme(4, :))
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('Coral cover (\% of reef area)', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title(['Total coral cover at reef ', num2str(reef)], ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('No control', '100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)
% 
% % Starfish at one reef ----------------------------------------------------
% figure(19), clf, hold on, grid on
% plot(t_vec_s0, N_y_2_s0(reef, :), 'Linewidth', 2, 'Color', 'black')
% plot(t_vec_s1, N_y_2_s1(reef, :), 'Linewidth', 2, 'Color', colour_scheme(1, :))
% plot(t_vec_s2, N_y_2_s2(reef, :), 'Linewidth', 2, 'Color', colour_scheme(2, :))
% plot(t_vec_s3, N_y_2_s3(reef, :), 'Linewidth', 2, 'Color', colour_scheme(3, :))
% plot(t_vec_s4, N_y_2_s4(reef, :), 'Linewidth', 2, 'Color', colour_scheme(4, :))
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 2+ starfish', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title(['No. of age 2+ starfish at reef ', num2str(reef)], ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('No control', '100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)
% 
% % Age 2+ starfish at one reef (control scenarios only) --------------------
% figure(20), clf, hold on, grid on
% plot(t_vec_s1, N_y_2_s1(reef, :), 'Linewidth', 2, 'Color', colour_scheme(1, :))
% % plot(t_vec_s2, N_y_2_s2(reef, :), 'Linewidth', 2, 'Color', colour_scheme(2, :))
% % plot(t_vec_s3, N_y_2_s3(reef, :), '--', 'Linewidth', 2, 'Color', colour_scheme(31, :))
% % plot(t_vec_s4, N_y_2_s4(reef, :), 'Linewidth', 2, 'Color', colour_scheme(4, :))
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 2+ starfish', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title(['No. of age 2+ starfish at reef ', num2str(reef)], ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% % legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
% %     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
% %     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)
% 
% % Age 1 starfish at one reef (control scenarios only) ---------------------
% figure(21), clf, hold on, grid on
% plot(t_vec_s1, N_y_1_s1(reef, :), 'Linewidth', 2, 'Color', colour_scheme(1, :))
% % plot(t_vec_s2, N_y_1_s2(reef, :), 'Linewidth', 2, 'Color', colour_scheme(2, :))
% % plot(t_vec_s3, N_y_1_s3(reef, :), '--', 'Linewidth', 2, 'Color', colour_scheme(3, :))
% % plot(t_vec_s4, N_y_1_s4(reef, :), 'Linewidth', 2, 'Color', colour_scheme(4, :))
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 1 starfish', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title(['No. of age 1 starfish at reef ', num2str(reef)], ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% % legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
% %     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
% %     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)
% 
% % Age 0 starfish at one reef (control scenarios only) ---------------------
% figure(22), clf, hold on, grid on
% plot(t_vec_s1, N_y_0_s1(reef, :), 'Linewidth', 2, 'Color', colour_scheme(1, :))
% % plot(t_vec_s2, N_y_0_s2(reef, :), 'Linewidth', 2, 'Color', colour_scheme(2, :))
% % plot(t_vec_s3, N_y_0_s3(reef, :), '--', 'Linewidth', 2, 'Color', colour_scheme(3, :))
% % plot(t_vec_s4, N_y_0_s4(reef, :), 'Linewidth', 2, 'Color', colour_scheme(4, :))
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 0 starfish', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title(['No. of age 0 starfish at reef ', num2str(reef)], ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% % legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
% %     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
% %     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)



%% STARFISH MAP ===========================================================
% % Starfish population on GBR (categories) ---------------------------------
% figure(23), clf, hold on, box on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% xlim([140, 155])
% ylim([-26, -8])
% % Plot reef locations by colour based on starfish presence
% for i = 1:length(lon)
%     if N_y_2_s4(i, end) >= 10^6
%         p7 = plot(lon(i), lat(i), '.', 'Markersize', 12, 'Color', magma_palette(1, :));
%     elseif N_y_2_s4(i, end) >= 10^4 && N_y_2_s4(i, end) < 10^6
%         p8 = plot(lon(i), lat(i), '.', 'Markersize', 12, 'Color', magma_palette(2, :));
%     elseif N_y_2_s4(i, end) >= 10^2 && N_y_2_s4(i, end) < 10^4
%         p9 = plot(lon(i), lat(i), '.', 'Markersize', 12, 'Color', magma_palette(3, :));
%     elseif N_y_2_s4(i, end) >= 0.5 && N_y_2_s4(i, end) < 10^2
%         p10 = plot(lon(i), lat(i), '.', 'Markersize', 12, 'Color', magma_palette(4, :));
%     else
%         p11 = plot(lon(i), lat(i), '.', 'Markersize', 12, 'Color', magma_palette(5, :));
%     end
% end
% % Add labels
% set(gca, 'FontSize', ticks_FS, 'BoxStyle', 'Full');
% title(['Adult starfish population after ', num2str(t_end), ' years with no control'], ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% [h, icons] = legend([p7 p8 p9 p10 p11], '$>$ 10$^6$ starfish', '10$^4-$10$^6$ starfish', '10$^2-$10$^4$ starfish', ...
%     '$<$ 10$^2$ starfish', '0 starfish', 'Interpreter', 'Latex', 'Fontsize', legend_FS);
% icons = findobj(icons, 'Type', 'line');
% icons = findobj(icons, 'Marker', 'none', '-xor');
% set(icons, 'MarkerSize', 25)

% % Starfish population on GBR (colormap) -----------------------------------
% % for i = 1:t_end+1
% figure(24), clf, hold on, box on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% % Plot reef locations by color depending on initial cots numbers
% scatter(lon, lat, 12, N_y_2_s0(:, t_end), 'filled')
% c2 = colorbar;
% colormap(magma)
% % caxis([0 1850])
% % Focus the figure on GBR and QLD
% xlim([140, 155])
% ylim([-26, -8])
% % Add labels
% set(gca, 'FontSize', ticks_FS, 'BoxStyle', 'Full');
% title(['Adult starfish population ', num2str(t_end), ' years after outbreak'], ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% c2.Label.String = 'No. of age 2+ starfish';
% c2.Label.Interpreter = 'Latex';
% c2.Label.FontSize = 14;
% % end



%% STARFISH SPREAD LATITUDE ===============================================

% Get adult starfish at final timestep ------------------------------------
starfish_final_s0 = N_y_2_s0(:, end);
% starfish_s0_per_sqkm = (starfish_final_s0./reef_area);

% Sort final starfish by latitude of reefs 
[lat_sorted, sort_index] = sort(lat);
starfish_s0_final_sorted = starfish_final_s0(sort_index);

% Initialise indexes for new arrays
index_500 = 1;

% Select data with reefs that have at least 500 starfish
for i = 1:num_reefs
    if starfish_s0_final_sorted(i) >= 500
        starfish_final_reefs(index_500) = starfish_s0_final_sorted(i);
        lat_starfish_values(index_500) = lat_sorted(i);
        index_500 = index_500 + 1;
    end
end

% Plot of starfish at the end of outbreak ---------------------------------
figure(25), clf

% Adult starfish on GBR (colormap) ----------------------------------------
p1 = subplot(1, 2, 1); hold on, box on
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% Plot reef locations by color depending on initial cots numbers
scatter(lon, lat, 12, starfish_final_s0, 'filled')
c2 = colorbar;
cmap = colormap(magma);
cmap = flipud(cmap);
colormap(cmap)
% Focus the figure on GBR and QLD
xlim([140, 155])
ylim([-25, -10])
p1.Position = [0.07 0.08 0.37 0.84];
% Add labels
set(gca, 'FontSize', ticks_FS, 'BoxStyle', 'Full');
xlabel('Longitude', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Latitude', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title(['Adult starfish population ', num2str(t_end), ' years after outbreak'], ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
c2.Label.String = 'No. of age 2+ starfish';
c2.Label.Interpreter = 'Latex';
c2.Label.FontSize = 14;

% Histogram of reefs with at least 1000 starfish --------------------------
p2 = subplot(1, 2, 2); hold on, box on
hs = histogram(lat_starfish_values, 40, 'Orientation', 'horizontal');
hs.BinWidth = 0.35;
set(gca, 'FontSize', ticks_FS)
p2.Position = [0.6 0.08 0.35 0.84];
ylim([-25, -10])
xlabel('No. of reefs with $\geq500$ adult starfish', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Latitude', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title({'Latitudinal spread of reefs with at least 500 ', ...
    ['adult starfish after ', num2str(t_end), ' years']}, ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)

% % Save - to avoid font resizing in colorbar
% saveas(gcf, 'Plots/03_paper/starfish_outbreak_map_histogram.png')



%% CORAL COVER MAP ========================================================

% % Coral cover on GBR (colormap) -------------------------------------------
% figure(26), clf, hold on, box on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% % Plot reef locations by color depending on initial cots numbers
% scatter(lon, lat, 12, C_y_f_s1(:, end), 'filled')
% % Focus the figure on GBR and QLD
% xlim([140, 155])
% ylim([-26, -9])
% % Add labels
% set(gca, 'FontSize', ticks_FS, 'BoxStyle', 'Full');
% title(['Coral cover after ', num2str(t_end), ' years with 100\% effort over 1737.1 $km^2$'], ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% % Add colorbar
% c = colorbar;
% cmap = colormap(viridis);
% cmap = flipud(cmap);
% colormap(cmap)
% % caxis([0 1])
% c.Label.String = 'Coral cover (\% of reef area)';
% c.Label.Interpreter = 'Latex';
% c.Label.FontSize = 14;
% 
% % % Save - to avoid font resizing in colorbar
% % saveas(gcf, 'Plots/04_comparisons/coral_gbr_map_compare.png')



%% CORAL HISTOGRAMS =======================================================
% % Histograms of coral cover -----------------------------------------------
% figure(27), clf, hold on
% for i = 1:5
%     subplot(2, 3, i), hold on
%     histogram(coral_end_norm(:, i))
% %     ylim([0 1400])
%     if i == 1
%         title('No control')
%     else
%         title(['Control scenario ', num2str(i)])
%     end
% end

% % Histograms of extra coral cover -----------------------------------------
% figure(28), clf, hold on
% for i = 1:size(coral_end_diff_norm, 2)
%     subplot(2, 2, i), hold on
%     histogram(coral_end_diff_norm(:, i))
% %     xlim([-0.05 0.55])
% %     ylim([0 2000])
%     title(['Control scenario ', num2str(i)])
% end



%% CORAL BAR GRAPH CATEGORICAL ============================================
% Bar graph of extra coral cover in categories ----------------------------
scenarios = categorical({'100\% over 1737.1 $km^2$', '50\% over 3474.2 $km^2$', ...
    '25\% over 6948.4 $km^2$', '12.2\% over 14240 $km^2$'});
scenarios = reordercats(scenarios, {'100\% over 1737.1 $km^2$', '50\% over 3474.2 $km^2$', ...
    '25\% over 6948.4 $km^2$', '12.2\% over 14240 $km^2$'});
figure(29), clf, hold on
b = bar(scenarios, coral_compare, 'EdgeColor', 'None');
b(1).FaceColor = viridis_palette_5(1, :);
b(2).FaceColor = viridis_palette_5(2, :);
b(3).FaceColor = viridis_palette_5(3, :);
b(4).FaceColor = viridis_palette_5(4, :);
b(5).FaceColor = viridis_palette_5(5, :);
set(gca, 'FontSize', ticks_FS+1);
set(gca, 'TickLabelInterpreter', 'Latex')
% ylim([-10 16])
title('\qquad\qquad\qquad Number of additional reefs with $x\%$ coral cover compared to no control after 100 years', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS);
xlabel('Control scenario', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
xh = get(gca,'xlabel');
px = get(xh,'position');
set(xh,'position',px)
ylabel('No. of additional reefs', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
legend('$x\leq 10$\% coral cover', '$10\%< x\leq30\%$ coral cover', ...
    '$30\%<x\leq 50\%$ coral cover', '$50\%<x\leq 75\%$ coral cover', ...
    '$x> 75\%$ coral cover', 'Interpreter', 'Latex', 'Fontsize', ...
    legend_FS, 'Location', 'NorthEastOutside');



%% CORAL COVER LATITUDE ===================================================
% % Histogram of reefs with at least 1% increase ----------------------------
% figure(30), clf, hold on
% h1 = histogram(lat_sorted, 'BinWidth', 0.3);
% h2 = histogram(lat_more_1_percent, 40);
% set(gca, 'XDir', 'Reverse')
% set(gca, 'FontSize', ticks_FS)
% xlabel('Latitude', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of reefs with $\geq1\%$ coral cover increase', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title({'Latitudinal spread of reefs with at least 1\% coral cover increase after 100 years', ...
%     'with 25\% control effort at 672 reefs compared to no control'}, ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)

% % Vertical spread of coral cover increase ---------------------------------
% figure(31), clf, hold on
% % lat(i) > -17 && lat(i) < -14.75
% patch('XData', [-17 -17 -14.75 -14.75], 'YData', [0.01 6 6 0.01], ...
%     'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'None')
% plot(lat_sorted, coral_best_compare_sorted, 'Linewidth', 1, 'Color', colour_scheme(3, :))
% ylim([0 6])
% xlim([min(lat_sorted), max(lat_sorted)])
% set(gca, 'XDir', 'Reverse')
% set(gca, 'FontSize', ticks_FS)
% xlabel('Latitude', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('Additional coral cover (\% of reef area)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title(['Additional coral cover after ', num2str(t_end), ' years with 25$\%$ effort at 672 reefs'], ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)



%% CORAL COVER INCREASE ===================================================

% Comparing percentage of coral cover -------------------------------------
% Calculate difference in coral cover percentage to no control 
coral_final_s0_percent = C_y_f_s0(:, end) ./ reef_area;
coral_final_s1_percent = C_y_f_s1(:, end) ./ reef_area;
coral_final_s2_percent = C_y_f_s2(:, end) ./ reef_area;
coral_final_s3_percent = C_y_f_s3(:, end) ./ reef_area;
coral_final_s4_percent = C_y_f_s4(:, end) ./ reef_area;
coral_s1_compare = (coral_final_s1_percent - coral_final_s0_percent) * 100;
coral_s2_compare = (coral_final_s2_percent - coral_final_s0_percent) * 100;
coral_s3_compare = (coral_final_s3_percent - coral_final_s0_percent) * 100;
coral_s4_compare = (coral_final_s4_percent - coral_final_s0_percent) * 100;

% Sort percentage increase by latitude of reefs 
[lat_sorted, sort_index] = sort(lat);
coral_s1_compare_sorted = coral_s1_compare(sort_index);
coral_s2_compare_sorted = coral_s2_compare(sort_index);
coral_s3_compare_sorted = coral_s3_compare(sort_index);
coral_s4_compare_sorted = coral_s4_compare(sort_index);

% Initialise indexes for new arrays
index_s1 = 1;
index_s2 = 1;
index_s3 = 1;
index_s4 = 1;

% Select data with reefs that have at least 1% coral cover increase
for i = 1:num_reefs
    % Compare scenario 1
    if coral_s1_compare_sorted(i) >= 1
        coral_s1_compare_1_percent(index_s1) = coral_s1_compare_sorted(i);
        lat_s1_1_percent(index_s1) = lat_sorted(i);
        index_s1 = index_s1 + 1;
    end
    % Compare scenario 2
    if coral_s2_compare_sorted(i) >= 1
        coral_s2_compare_1_percent(index_s2) = coral_s2_compare_sorted(i);
        lat_s2_1_percent(index_s2) = lat_sorted(i);
        index_s2 = index_s2 + 1;
    end
    % Compare scenario 3
    if coral_s3_compare_sorted(i) >= 1
        coral_s3_compare_1_percent(index_s3) = coral_s3_compare_sorted(i);
        lat_s3_1_percent(index_s3) = lat_sorted(i);
        index_s3 = index_s3 + 1;
    end
    % Compare scenario 4
    if coral_s4_compare_sorted(i) >= 1
        coral_s4_compare_1_percent(index_s4) = coral_s4_compare_sorted(i);
        lat_s4_1_percent(index_s4) = lat_sorted(i);
        index_s4 = index_s4 + 1;
    end
end

% Number of reefs with increase 1% increase
index_s1
index_s2
index_s3
index_s4

% Plot of coral cover increase --------------------------------------------
figure(32), clf

% Coral cover on GBR (colormap) -------------------------------------------
sp1 = subplot(1, 2, 1); hold on, box on
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% Plot reef locations by color depending on initial cots numbers
scatter(lon, lat, 12, coral_s1_compare, 'filled')
% Focus the figure on GBR and QLD
xlim([140, 155])
ylim([-25, -10])
sp1.Position = [0.07 0.08 0.43 0.84];
% Add labels
set(gca, 'FontSize', ticks_FS, 'BoxStyle', 'Full');
xlabel('Longitude', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Latitude', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title({['Coral cover increase after ', num2str(t_end), ' years with 100\% effort'], ...
    'over 1737.1 $km^2$ compared to no control'}, ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
% Add colorbar
c = colorbar;
cmap = colormap(viridis);
cmap = flipud(cmap);
colormap(cmap)
c.Label.String = 'Coral cover increase (\% of reef area)';
c.Label.Interpreter = 'Latex';
c.Label.FontSize = 14;

% Histogram for all four scenarios ----------------------------------------
sp2 = subplot(1, 2, 2); hold on

% Plot histogram first and get info
h1 = histogram(lat_s1_1_percent, 40, 'Orientation', 'horizontal');
h2 = histogram(lat_s2_1_percent, 40, 'Orientation', 'horizontal');
h3 = histogram(lat_s3_1_percent, 40, 'Orientation', 'horizontal');
h4 = histogram(lat_s4_1_percent, 40, 'Orientation', 'horizontal');
h1.BinWidth = 0.5;
h2.BinWidth = 0.5;
h3.BinWidth = 0.5;
h4.BinWidth = 0.5;

% Get the values to then make a scatter plot with
% x axis values - number of reefs
s1_reefs = h1.Values;
s2_reefs = h2.Values;
s3_reefs = h3.Values;
s4_reefs = h4.Values;

% y axis values
lat_values_s1 = h1.BinEdges + (h1.BinWidth/2);
lat_values_s1 = lat_values_s1(2:end);
lat_values_s2 = h2.BinEdges + (h2.BinWidth/2);
lat_values_s2 = lat_values_s2(2:end);
lat_values_s3 = h3.BinEdges + (h3.BinWidth/2);
lat_values_s3 = lat_values_s3(2:end);
lat_values_s4 = h4.BinEdges + (h4.BinWidth/2);
lat_values_s4 = lat_values_s4(2:end);

% Create arrays and delete rows with zeros
s1_reefs_plot = [s1_reefs; lat_values_s1]';
s1_reefs_plot = s1_reefs_plot(all(s1_reefs_plot, 2), :);
s2_reefs_plot = [s2_reefs; lat_values_s2]';
s2_reefs_plot = s2_reefs_plot(all(s2_reefs_plot, 2), :);
s3_reefs_plot = [s3_reefs; lat_values_s3]';
s3_reefs_plot = s3_reefs_plot(all(s3_reefs_plot, 2), :);
s4_reefs_plot = [s4_reefs; lat_values_s4]';
s4_reefs_plot = s4_reefs_plot(all(s4_reefs_plot, 2), :);

% Clear this subplot
cla(sp2)

% Now make the scatter plot
s1 = scatter(s1_reefs_plot(:, 1), s1_reefs_plot(:, 2), 'o');
s1.MarkerFaceColor = colour_scheme(1, :);
s1.MarkerEdgeColor = colour_scheme(1, :);
s1.SizeData = 70;
s1.LineWidth = 2;

s2 = scatter(s2_reefs_plot(:, 1), s2_reefs_plot(:, 2), 'x');
s2.MarkerFaceColor = colour_scheme(2, :);
s2.MarkerEdgeColor = colour_scheme(2, :);
s2.SizeData = 90;
s2.LineWidth = 2;

s3 = scatter(s3_reefs_plot(:, 1), s3_reefs_plot(:, 2), '+');
s3.MarkerFaceColor = colour_scheme(3, :);
s3.MarkerEdgeColor = colour_scheme(3, :);
s3.SizeData = 90;
s3.LineWidth = 2;

s4 = scatter(s4_reefs_plot(:, 1), s4_reefs_plot(:, 2), 'o');
s4.MarkerFaceColor = colour_scheme(4, :);
s4.MarkerEdgeColor = colour_scheme(4, :);
s4.SizeData = 70;
s4.LineWidth = 2;

% Set the axis, labels etc. 
set(gca, 'FontSize', ticks_FS)
sp2.Position = [0.6 0.08 0.35 0.84];
ylim([-25, -10])
% xlim([0 60])
xlabel('No. of reefs with $\geq 1\%$ coral cover increase', ...
    'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Latitude', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title({'Latitudinal spread of reefs with at least 1\%', ...
    ['coral cover increase after ', num2str(t_end), ' years']}, ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
legend('100\% effort over 1737.1 $km^2$', '50\% effort over 3474.2 $km^2$', ...
    '25\% effort over 6948.4 $km^2$', '12.2\% effort over 14240 $km^2$', ...
    'Interpreter', 'Latex', 'Fontsize', legend_FS, 'Location', 'SouthEast');

% % Save - to avoid font resizing in colorbar
% saveas(gcf, 'Plots/03_paper/coral_cover_increase_map_histogram_all.png')


% Plot of coral cover increase (best scenario only) -----------------------
figure(33), clf

% Coral cover on GBR (colormap) -------------------------------------------
sp3 = subplot(1, 2, 1); hold on, box on
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% Plot reef locations by color depending on initial cots numbers
scatter(lon, lat, 12, coral_s1_compare, 'filled')
% Focus the figure on GBR and QLD
xlim([140, 155])
ylim([-25, -10])
sp3.Position = [0.07 0.08 0.43 0.84];
% Add labels
set(gca, 'FontSize', ticks_FS, 'BoxStyle', 'Full');
xlabel('Longitude', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Latitude', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title({['Coral cover increase after ', num2str(t_end), ' years with 100\% effort'], ...
    'over 1737.1 $km^2$ compared to no control'}, ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
% Add colorbar
c = colorbar;
cmap = colormap(viridis);
cmap = flipud(cmap);
colormap(cmap)
c.Label.String = 'Coral cover increase (\% of reef area)';
c.Label.Interpreter = 'Latex';
c.Label.FontSize = 14;

% Histogram for best scenario only ----------------------------------------
sp4 = subplot(1, 2, 2); hold on
h5 = histogram(lat_s1_1_percent, 40, 'Orientation', 'horizontal');
h5.BinWidth = 0.35;
set(gca, 'FontSize', ticks_FS)
sp4.Position = [0.6 0.08 0.35 0.84];
ylim([-25, -10])
xlabel('No. of reefs with $\geq 1\%$ coral cover increase', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Latitude', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title({'Latitudinal spread of reefs with at least 1\%', ...
    ['coral cover increase after ', num2str(t_end), ' years']}, ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)


% % Save - to avoid font resizing in colorbar
% saveas(gcf, 'Plots/03_paper/coral_cover_increase_map_histogram_scenario1.png')


%% CONTROL EFFORT =========================================================
% Control effort heatmaps -------------------------------------------------
figure(34), clf, hold on, box on
label_strings = {'(a) 100\% effort over 1737.1 $km^2$', '(b) 50\% effort over 3474.2 $km^2$', ...
    '(c) 25\% effort over 6948.4 $km^2$', '(d) 12.2\% effort over 14240 $km^2$'};
for i = 1:size(control_plot, 2)
    sp = subplot(2, 2, i); hold on, box on
    % Plot outline of Australia
    pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
    % Plot reef locations by color depending on control effort
    scatter(lon, lat, 10, control_plot(:, i), 'filled')
    colorbar off
    colormap parula
    caxis([0 1])
    % Focus the figure on GBR and QLD
    xlim([140, 155])
    ylim([-25, -10])
    set(gca, 'YTick', -25:5:-10);
    % Add labels
    set(gca, 'BoxStyle', 'Full');
    xlabel(label_strings(i), 'Interpreter', 'Latex', 'Fontsize', axis_FS-2)
    title(['Control scenario ', num2str(i)], 'Interpreter', 'Latex', 'Fontsize', title_FS-3)
    % Positioning
    if i == 1
        sp.Position = [0.07 0.58 0.33 0.35];
    elseif i == 2
        sp.Position = [0.5 0.58 0.33 0.35];
    elseif i == 3
        sp.Position = [0.07 0.11 0.33 0.35];
    elseif i == 4
        sp.Position = [0.5 0.11 0.33 0.35];
    end
end

% Add colorbar
c = colorbar;
c.Position = [0.88 0.11 0.02 0.815];
c.Label.String = 'Percentage of adult starfish culled ($k_{i,t}$)';
c.Label.Interpreter = 'Latex';
c.Label.FontSize = 14;

% % Save - to avoid font resizing in colorbar
% saveas(gcf, 'Plots/03_paper/control_effort_compare.png')

% % Control effort heatmaps (simple) ----------------------------------------
% figure(34), clf, hold on, box on
% label_strings = {'(a) 100\% starfish culled at 168 reefs', '(b) 50\% starfish culled at 336 reefs', ...
%     '(c) 25\% starfish culled at 672 reefs', '(d) 7.72\% starfish culled at 2175 reefs'};
% for i = 1:size(control_plot, 2)
%     sp = subplot(2, 2, i); hold on, box on
%     % Plot outline of Australia
%     pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
%     % Plot reef locations by color depending on control effort
%     for j = 1:size(control_plot, 1)
%         if control_plot(j, i) > 0
%             p1 = plot(lon(j), lat(j), '.', 'MarkerSize', 12, 'Color', viridis_palette_4(2, :));
%         else
%             p2 = plot(lon(j), lat(j), '.', 'MarkerSize', 12, 'Color', viridis_palette_4(3, :));
%         end
%     end
%     % Focus the figure on GBR and QLD
%     xlim([140, 155])
%     ylim([-26, -8])
%     set(gca, 'YTick', -26:6:-8);
%     % Add labels
%     set(gca, 'BoxStyle', 'Full');
%     xlabel(label_strings(i), 'Interpreter', 'Latex', 'Fontsize', axis_FS-2)
%     title(['Strategy ', num2str(i)], 'Interpreter', 'Latex', 'Fontsize', title_FS-3)
%     % Legend for first plot only
%     if i == 1
%         [h, icons] = legend([p1 p2], 'Starfish culled', 'No starfish culled', 'Interpreter', 'Latex', 'Fontsize', title_FS-2);
%         icons = findobj(icons, 'Type', 'line');
%         icons = findobj(icons, 'Marker', 'none', '-xor');
%         set(icons, 'MarkerSize', 20)
%     end
% end
