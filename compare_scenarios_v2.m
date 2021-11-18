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


%% SCENARIO 0: No control, but simulation run with same conditions

% CONTROL EFFORT ----------------------------------------------------------
% No control
control_effort_s0 = 0;

% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s0, C_y_f_s0, N_y_2_s0, N_y_1_s0, N_y_0_s0, tau_ratio_s0] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s0, dispersal_eq);

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
% Cull all the starfish only in the initiation box
control_effort_s1 = zeros(num_reefs, t_end);
for i = 1:num_reefs
    if (lat(i) > -17 && lat(i) < -14.75) && (lon(i) > 145 && lon(i) < 147)
        control_effort_s1(i, :) = 1;
    end
end

% Count number of reefs controlled in control scenario 1
num_reefs_control_s1 = nnz(control_effort_s1(:, 1))
budget_s1 = sum(control_effort_s1, 'all')

% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s1, C_y_f_s1, N_y_2_s1, N_y_1_s1, N_y_0_s1, tau_ratio_s1] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s1, dispersal_eq);

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
    if (lat(i) > -18.45 && lat(i) < -14) && (lon(i) > 143 && lon(i) < 147.75)
        control_effort_s2(i, :) = 0.5;
    end
end

% Print the total number of reefs and 'budget' for control scenario
num_reefs_control_s2 = nnz(control_effort_s2(:, 1))
budget_s2 = sum(control_effort_s2, 'all')

% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s2, C_y_f_s2, N_y_2_s2, N_y_1_s2, N_y_0_s2, tau_ratio_s2] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s2, dispersal_eq);

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
    if (lat(i) > -20 && lat(i) < -12.69) && (lon(i) > 142 && lon(i) < 152)
        control_effort_s3(i, :) = 0.25;
    end
end

% Print the total number of reefs and 'budget' for control scenario
num_reefs_control_s3 = nnz(control_effort_s3(:, 1))
budget_s3 = sum(control_effort_s3, 'all')

% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s3, C_y_f_s3, N_y_2_s3, N_y_1_s3, N_y_0_s3, tau_ratio_s3] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s3, dispersal_eq);

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
% Cull starfish everywhere based on total budget
effort_per_year = budget_s1 / t_end;
effort_per_reef = effort_per_year / num_reefs
control_effort_s4 = effort_per_reef * ones(num_reefs, t_end);

% Print the total number of reefs and 'budget' for control scenario
num_reefs_control_s4 = nnz(control_effort_s4(:, 1))
budget_s4 = sum(control_effort_s4, 'all')

% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_s4, C_y_f_s4, N_y_2_s4, N_y_1_s4, N_y_0_s4, tau_ratio_s4] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s4, dispersal_eq);

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
% Cull all the starfish only in the initiation box for all 50 years
control_effort_other = ones(num_reefs, t_end);

% Count number of reefs controlled in control scenario 1
num_reefs_control_s1 = nnz(control_effort_other(:, 1))
budget_other = sum(control_effort_other, 'all')

% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec_other, C_y_f_other, N_y_2_other, N_y_1_other, N_y_0_other, tau_ratio_other] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_other, dispersal_eq);

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
coral_end_s0 = C_y_f_s0(:, end);
coral_end_s1 = C_y_f_s1(:, end);
coral_end_s2 = C_y_f_s2(:, end);
coral_end_s3 = C_y_f_s3(:, end);
coral_end_s4 = C_y_f_s4(:, end);

% Combine into matrix
coral_end = [coral_end_s0 coral_end_s1 coral_end_s2 coral_end_s3 coral_end_s4];

% Normalise
coral_end_norm = normalize(coral_end, 'center');

% Initialise matrix to store comparisons
coral_compare_all = zeros(size(coral_end, 2), 4);
% coral_compare = zeros(size(coral_end, 2)-1, 4);

% Loop over each column in matrix and count number of reefs
for i = 1:size(coral_end, 2)
    % Count the number of reefs with less than 1% coral 
    coral_compare_all(i, 1) = sum(coral_end(:, i) < 0.01);
    
    % Count the number of reefs with between 1% and 5% coral
    coral_compare_all(i, 2) = sum(coral_end(:, i) >= 0.01 & coral_end(:, i) < 0.05);
    
    % Count the number of reefs with between 5% and 30% coral
    coral_compare_all(i, 3) = sum(coral_end(:, i) >= 0.05 & coral_end(:, i) < 0.3);
    
    % Count the number of reefs with more than 30%
    coral_compare_all(i, 4) = sum(coral_end(:, i) >= 0.3);
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


%% PLOTS
% Fontsizes for plotting
axis_FS = 15;
title_FS = 16;
legend_FS = 13;
ticks_FS = 12;

% Colours for plotting
colour_scheme = cbrewer('div', 'RdBu', 4);
colour_scheme_brown = cbrewer('div', 'BrBG', 9);
colour_scheme_purple = cbrewer('div', 'PRGn', 9);

% Define viridis palette colours 
viridis_palette_4 = [253, 231, 37; 53, 183, 121; 49, 104, 142; 68, 1, 84];
viridis_palette_4 = viridis_palette_4/255;


%% CORAL HISTOGRAMS =======================================================
% Histogram of coral cover ------------------------------------------------
% figure(5), clf, hold on
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

% % Histogram of extra coral cover ------------------------------------------
% figure(10), clf, hold on
% for i = 1:size(coral_end_diff_norm, 2)
%     subplot(2, 2, i), hold on
%     histogram(coral_end_diff_norm(:, i))
% %     xlim([-0.05 0.55])
% %     ylim([0 2000])
%     title(['Control scenario ', num2str(i)])
% end

% Bar graph of extra coral cover in categories ----------------------------
scenarios = categorical({'100\% at 168 reefs', '50\% at 336 reefs', ...
    '25\% at 672 reefs', '7.72\% at 2175 reefs'});
scenarios = reordercats(scenarios, {'100\% at 168 reefs', '50\% at 336 reefs', ...
    '25\% at 672 reefs', '7.72\% at 2175 reefs'});
figure(5), clf, hold on
b = bar(scenarios, coral_compare, 'EdgeColor', 'None');
b(1).FaceColor = viridis_palette_4(1, :);
b(2).FaceColor = viridis_palette_4(2, :);
b(3).FaceColor = viridis_palette_4(3, :);
b(4).FaceColor = viridis_palette_4(4, :);
set(gca, 'FontSize', ticks_FS+1);
set(gca, 'TickLabelInterpreter', 'Latex')
title('\qquad\qquad\qquad Number of additional reefs with $x\%$ coral cover compared to no control after 100 years', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS);
xlabel('Control scenario', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
xh = get(gca,'xlabel');
px = get(xh,'position');
px(2) = -1 + px(2);
set(xh,'position',px)
ylabel('No. of additional reefs', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
legend('$x<1$\% coral cover', '$1\%\leq x<5\%$ coral cover', ...
    '$5\%\leq x<30\%$ coral cover', '$x\geq 30\%$ coral cover', ...
    'Interpreter', 'Latex', 'Fontsize', legend_FS, 'Location', 'NorthEastOutside');


%% CORAL ==================================================================
% GBR total coral cover ---------------------------------------------------
figure(1), clf, hold on
yline(num_reefs, '--', 'Linewidth', 2, 'Color', [0.5 0.5 0.5])
plot(t_vec_s0, coral_s0, 'Linewidth', 2, 'Color', 'black')
plot(t_vec_other, coral_other, '-.', 'Linewidth', 2, 'Color', colour_scheme_brown(1, :))
plot(t_vec_s1, coral_s1, 'Linewidth', 2, 'Color', colour_scheme(1, :))
plot(t_vec_s2, coral_s2, 'Linewidth', 2, 'Color', colour_scheme(2, :))
plot(t_vec_s3, coral_s3, 'Linewidth', 2, 'Color', colour_scheme(3, :))
plot(t_vec_s4, coral_s4, 'Linewidth', 2, 'Color', colour_scheme(4, :))
set(gca, 'FontSize', ticks_FS);
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total coral cover', 'Interpreter', 'Latex', ...
    'Fontsize', axis_FS)
title('Total coral cover on GBR', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
legend('Coral carrying capacity', 'No control', '100\% effort at 2175 reefs', ...
    '100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
    '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs',  ...
    'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% GBR total coral cover (additional) --------------------------------------
figure(2), clf, hold on
plot(t_vec_s1, coral_s1-coral_s0, 'Linewidth', 2, 'Color', colour_scheme(1, :))
plot(t_vec_s2, coral_s2-coral_s0, 'Linewidth', 2, 'Color', colour_scheme(2, :))
plot(t_vec_s3, coral_s3-coral_s0, 'Linewidth', 2, 'Color', colour_scheme(3, :))
plot(t_vec_s4, coral_s4-coral_s0, 'Linewidth', 2, 'Color', colour_scheme(4, :))
set(gca, 'FontSize', ticks_FS);
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total coral cover increase', 'Interpreter', 'Latex', ...
    'Fontsize', axis_FS)
title('Total coral cover increase on GBR compared to no control', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
    '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
    'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Initiation box total coral cover ----------------------------------------
% figure(3), clf, hold on
% yline(nnz(initial_state.N_0_2), '--', 'Linewidth', 2, 'Color', [0.5 0.5 0.5])
% plot(t_vec_s0, coral_box_s0, 'Linewidth', 2, 'Color', 'black')
% plot(t_vec_s1, coral_box_s1, 'Linewidth', 2, 'Color', colour_scheme(1, :))
% plot(t_vec_s2, coral_box_s2, 'Linewidth', 2, 'Color', colour_scheme(2, :))
% plot(t_vec_s3, coral_box_s3, 'Linewidth', 2, 'Color', colour_scheme(3, :))
% plot(t_vec_s4, coral_box_s4, 'Linewidth', 2, 'Color', colour_scheme(4, :))
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('Total coral cover', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title('Total coral cover in initiation box', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('Coral carrying capacity', 'No control', '100\% effort at 168 reefs', ...
%     '50\% effort at 336 reefs', '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% Initiation box total coral cover (additional) ---------------------------
figure(4), clf, hold on
plot(t_vec_s1, coral_box_s1-coral_box_s0, 'Linewidth', 2, 'Color', colour_scheme(1, :))
plot(t_vec_s2, coral_box_s2-coral_box_s0, 'Linewidth', 2, 'Color', colour_scheme(2, :))
plot(t_vec_s3, coral_box_s3-coral_box_s0, 'Linewidth', 2, 'Color', colour_scheme(3, :))
plot(t_vec_s4, coral_box_s4-coral_box_s0, 'Linewidth', 2, 'Color', colour_scheme(4, :))
set(gca, 'FontSize', ticks_FS);
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total coral cover increase', 'Interpreter', 'Latex', ...
    'Fontsize', axis_FS)
title('\qquad\qquad Total coral cover increase at initiation box compared to no control', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
    '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
    'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)


%% STARFISH ===============================================================
% GBR total age 2+ starfish -----------------------------------------------
% figure(5), clf, hold on
% plot(t_vec_s0, starfish_age2_s0, 'Linewidth', 2)
% plot(t_vec_s1, starfish_age2_s1, 'Linewidth', 2)
% plot(t_vec_s2, starfish_age2_s2, 'Linewidth', 2)
% plot(t_vec_s3, starfish_age2_s3, 'Linewidth', 2)
% plot(t_vec_s4, starfish_age2_s4, 'Linewidth', 2)
% set(gca, 'FontSize', ticks_FS)
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of adult starfish', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title('Total adult starfish population on GBR', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('No control', '100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthWest', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% GBR total age 2+ starfish (log scale) -----------------------------------
figure(6), clf, hold on
plot(t_vec_s0, starfish_age2_s0, 'Linewidth', 2, 'Color', 'Black')
plot(t_vec_s1, starfish_age2_s1, 'Linewidth', 2, 'Color', colour_scheme(1, :))
plot(t_vec_s2, starfish_age2_s2, 'Linewidth', 2, 'Color', colour_scheme(2, :))
plot(t_vec_s3, starfish_age2_s3, 'Linewidth', 2, 'Color', colour_scheme(3, :))
plot(t_vec_s4, starfish_age2_s4, 'Linewidth', 2, 'Color', colour_scheme(4, :))
set(gca, 'YScale', 'log')
set(gca, 'FontSize', ticks_FS)
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total no. of adult starfish (log scale)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title('Total adult starfish population on GBR', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
legend('No control', '100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
    '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
    'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% GBR decrease in age 2+ starfish (log scale) -----------------------------
figure(7), clf, hold on
plot(t_vec_s1, starfish_age2_s0-starfish_age2_s1, 'Linewidth', 2, 'Color', colour_scheme(1, :))
plot(t_vec_s2, starfish_age2_s0-starfish_age2_s2, 'Linewidth', 2, 'Color', colour_scheme(2, :))
plot(t_vec_s3, starfish_age2_s0-starfish_age2_s3, 'Linewidth', 2, 'Color', colour_scheme(3, :))
plot(t_vec_s4, starfish_age2_s0-starfish_age2_s4, 'Linewidth', 2, 'Color', colour_scheme(4, :))
set(gca, 'YScale', 'log')
set(gca, 'FontSize', ticks_FS)
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total adult starfish decrease (log scale)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title('\qquad\qquad\qquad\qquad\qquad\qquad Total decrease in adult starfish population on GBR compared to no control', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS-1)
legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
    '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
    'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % GBR total age 2+ starfish (control scenarios only) --------------------
% figure(7), clf, hold on, grid on
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

% % Initiation box total age 2+ starfish ------------------------------------
% figure(8), clf, hold on
% plot(t_vec_s0, starfish_age2_box_s0, 'Linewidth', 2)
% plot(t_vec_s1, starfish_age2_box_s1, 'Linewidth', 2)
% plot(t_vec_s2, starfish_age2_box_s2, 'Linewidth', 2)
% plot(t_vec_s3, starfish_age2_box_s3, 'Linewidth', 2)
% plot(t_vec_s4, starfish_age2_box_s4, 'Linewidth', 2)
% set(gca, 'FontSize', ticks_FS)
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of adult starfish', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title('Total adult starfish population at initiation box', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('No control', '100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthWest', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Initiation box total age 2+ starfish (log scale) ------------------------
% figure(8), clf, hold on
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

% Initiation box decrease in age 2+ starfish (log scale) ------------------
figure(9), clf, hold on
plot(t_vec_s1, starfish_age2_box_s0-starfish_age2_box_s1, 'Linewidth', 2, 'Color', colour_scheme(1, :))
plot(t_vec_s2, starfish_age2_box_s0-starfish_age2_box_s2, 'Linewidth', 2, 'Color', colour_scheme(2, :))
plot(t_vec_s3, starfish_age2_box_s0-starfish_age2_box_s3, 'Linewidth', 2, 'Color', colour_scheme(3, :))
plot(t_vec_s4, starfish_age2_box_s0-starfish_age2_box_s4, 'Linewidth', 2, 'Color', colour_scheme(4, :))
set(gca, 'YScale', 'log')
set(gca, 'FontSize', ticks_FS)
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total adult starfish decrease (log scale)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title('\quad\qquad\qquad\qquad\qquad\qquad\qquad\qquad Total decrease in adult starfish population at initiation box compared to no control', ...
    'Interpreter', 'Latex', 'Fontsize', title_FS-1)
legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
    '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
    'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Initiation box total age 2+ starfish (additional) -----------------------
% figure(10), clf, hold on, grid on
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
% figure(11), clf, hold on, grid on
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
% figure(12), clf, hold on, grid on
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
% figure(13), clf, hold on, grid on
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

% % Initiation box total starfish (control scenarios only) ------------
% figure(14), clf, hold on, grid on
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
% figure(15), clf, hold on, grid on
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


%% SINGLE REEF ============================================================
% % Coral cover at one reef -------------------------------------------------
% figure(16), clf, hold on, grid on
% plot(t_vec_s0, C_y_f_s0(1000, :), 'Linewidth', 2)
% plot(t_vec_s1, C_y_f_s1(1000, :), 'Linewidth', 2)
% plot(t_vec_s2, C_y_f_s2(1000, :), 'Linewidth', 2)
% plot(t_vec_s3, C_y_f_s3(1000, :), 'Linewidth', 2)
% plot(t_vec_s4, C_y_f_s4(1000, :), 'Linewidth', 2)
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('Coral cover (\% of reef area)', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title('Total coral cover at reef 1000', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('No control', '100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Starfish at one reef ----------------------------------------------------
% figure(17), clf, hold on, grid on
% plot(t_vec_s0, N_y_2_s0(1000, :), 'Linewidth', 2)
% plot(t_vec_s1, N_y_2_s1(1000, :), 'Linewidth', 2)
% plot(t_vec_s2, N_y_2_s2(1000, :), 'Linewidth', 2)
% plot(t_vec_s3, N_y_2_s3(1000, :), 'Linewidth', 2)
% plot(t_vec_s4, N_y_2_s4(1000, :), 'Linewidth', 2)
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 2+ starfish', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title('No. of age 2+ starfish at reef 1000', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend('No control', '100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
%     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
%     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Age 2+ starfish at one reef (control scenarios only) --------------------
% figure(18), clf, hold on, grid on
% plot(t_vec_s1, N_y_2_s1(1000, :), 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980])
% % plot(t_vec_s2, N_y_2_s2(1000, :), 'Linewidth', 2, 'Color', [0.9290 0.6940 0.1250])
% % plot(t_vec_s3, N_y_2_s3(1000, :), '--', 'Linewidth', 2, 'Color', [0.4940 0.1840 0.5560])
% % plot(t_vec_s4, N_y_2_s4(1000, :), 'Linewidth', 2, 'Color', [0.4660 0.6740 0.1880])
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 2+ starfish', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title('No. of age 2+ starfish at reef 1000', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% % legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
% %     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
% %     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Age 1 starfish at one reef (control scenarios only) ---------------------
% figure(19), clf, hold on, grid on
% plot(t_vec_s1, N_y_1_s1(1000, :), 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980])
% % plot(t_vec_s2, N_y_1_s2(1000, :), 'Linewidth', 2, 'Color', [0.9290 0.6940 0.1250])
% % plot(t_vec_s3, N_y_1_s3(1000, :), '--', 'Linewidth', 2, 'Color', [0.4940 0.1840 0.5560])
% % plot(t_vec_s4, N_y_1_s4(1000, :), 'Linewidth', 2, 'Color', [0.4660 0.6740 0.1880])
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 1 starfish', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title('No. of age 1 starfish at reef 1000', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% % legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
% %     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
% %     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Age 0 starfish at one reef (control scenarios only) ---------------------
% figure(20), clf, hold on, grid on
% plot(t_vec_s1, N_y_0_s1(1000, :), 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980])
% % plot(t_vec_s2, N_y_0_s2(1000, :), 'Linewidth', 2, 'Color', [0.9290 0.6940 0.1250])
% % plot(t_vec_s3, N_y_0_s3(1000, :), '--', 'Linewidth', 2, 'Color', [0.4940 0.1840 0.5560])
% % plot(t_vec_s4, N_y_0_s4(1000, :), 'Linewidth', 2, 'Color', [0.4660 0.6740 0.1880])
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 0 starfish', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title('No. of age 0 starfish at reef 1000', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% % legend('100\% effort at 168 reefs', '50\% effort at 336 reefs', ...
% %     '25\% effort at 672 reefs', '7.72\% effort at 2175 reefs', ...
% %     'Location', 'NorthEastOutside', 'Interpreter', 'Latex', 'Fontsize', legend_FS)

%% CORAL COVER MAP ========================================================
% Map calculations --------------------------------------------------------
coral_best_compare = (C_y_f_s4(:, end)-C_y_f_s0(:, end))*100;


% Coral cover on GBR (colormap) -------------------------------------------
figure(22), clf, hold on, box on
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% Plot reef locations by color depending on initial cots numbers
scatter(lon, lat, 12, coral_best_compare, 'filled')
% Focus the figure on GBR and QLD
xlim([140, 155])
ylim([-26, -8])
% Add labels
set(gca, 'FontSize', ticks_FS, 'BoxStyle', 'Full');
title(['Additional coral cover after ', num2str(t_end), ' years with 25$\%$ effort at 672 reefs'], ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
% Add colorbar
c = colorbar;
cmap = colormap(viridis);
cmap = flipud(cmap);
colormap(cmap)
c.Label.String = 'Additional coral cover (\% of reef area)';
c.Label.Interpreter = 'Latex';
c.Label.FontSize = 14;

% % Save - to avoid font resizing in colorbar
% saveas(gcf, 'Plots/04_comparisons/coral_gbr_map_compare.png')


%% CORAL COVER LATITUDE ===================================================
% Plot calculations -------------------------------------------------------
[lat_sorted, sort_index] = sort(lat);
coral_best_compare_sorted = coral_best_compare(sort_index);


% num_plot = 435;
% num_reef_per_point = num_reefs/num_plot;
% 
% lat_sorted_plot = lat_sorted(3:5:2173);
% coral_best_compare_sorted_plot = zeros(1, num_plot);
% 
% index = 1;
% 
% for i = 1:num_plot
%     this_sum = 0;
%     for j = 1:num_reef_per_point
%         this_sum = this_sum + coral_best_compare_sorted(index);
%         index = index + 1;
%     end
%     coral_best_compare_sorted_plot(i) = this_sum;
% end

% Select data with reefs that have at least 1\% coral cover increase
index = 1;
for i = 1:num_reefs
    if coral_best_compare_sorted(i) >= 1
        coral_best_compare_more(index) = coral_best_compare_sorted(i);
        lat_more(index) = lat_sorted(i);
        index = index + 1;
    end
end

% Histogram of reefs with at least 1% increase ----------------------------
figure(24), clf, hold on
% h1 = histogram(lat_sorted, 'BinWidth', 0.3);
h2 = histogram(lat_more, 40);
set(gca, 'XDir', 'Reverse')
set(gca, 'FontSize', ticks_FS)
xlabel('Latitude', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('No. of reefs with $\geq1\%$ coral cover increase', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title({'Latitudinal spread of reefs with at least 1\% coral cover increase after 100 years', ...
    'with 25\% control effort at 672 reefs compared to no control'}, ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)

% % Calculate density of reefs
% density_reefs = h2.Values ./ h1.Values(3:end);
% density_reefs(isnan(density_reefs)) = 0;

% % Plot density of reefs instead now ---------------------------------------
% figure(25), clf, hold on
% bar(density_reefs, 'BarWidth', 1)
% set(gca, 'XDir', 'Reverse')
% set(gca, 'FontSize', ticks_FS)
% xlabel('Latitude', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('Proportion of reefs with $\geq$1\% coral cover increase', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title('Latitudinal spread of reefs with $\geq$1\% coral cover increase', ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)


% % Vertical spread of coral cover increase ---------------------------------
% figure(23), clf, hold on
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


%% CONTROL EFFORT =========================================================
% % Control effort heatmaps -------------------------------------------------
% figure(21), clf, hold on, box on
% label_strings = {'(a) 100\% effort at 168 reefs', '(b) 50\% effort at 336 reefs', ...
%     '(c) 25\% effort at 672 reefs', '(d) 7.72\% effort at 2175 reefs'};
% for i = 1:size(control_plot, 2)
%     sp = subplot(2, 2, i); hold on, box on
%     % Plot outline of Australia
%     pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
%     % Plot reef locations by color depending on control effort
%     scatter(lon, lat, 10, control_plot(:, i), 'filled')
%     colorbar off
%     colormap parula
%     caxis([0 1])
%     % Focus the figure on GBR and QLD
%     xlim([140, 155])
%     ylim([-26, -8])
%     set(gca, 'YTick', -26:6:-8);
%     % Add labels
%     set(gca, 'BoxStyle', 'Full');
%     xlabel(label_strings(i), 'Interpreter', 'Latex', 'Fontsize', axis_FS-2)
%     title(['Control scenario ', num2str(i)], 'Interpreter', 'Latex', 'Fontsize', title_FS-3)
%     % Positioning
%     if i == 1
%         sp.Position = [0.07 0.58 0.33 0.35];
%     elseif i == 2
%         sp.Position = [0.5 0.58 0.33 0.35];
%     elseif i == 3
%         sp.Position = [0.07 0.11 0.33 0.35];
%     elseif i == 4
%         sp.Position = [0.5 0.11 0.33 0.35];
%     end
% end
% 
% % Add colorbar
% c = colorbar;
% c.Position = [0.88 0.11 0.02 0.815];
% c.Label.String = 'Percentage of adult starfish culled';
% c.Label.Interpreter = 'Latex';
% c.Label.FontSize = 14;
% 
% % % Save - to avoid font resizing in colorbar
% % saveas(gcf, 'Plots/04_comparisons/control_effort_compare.png')

figure(22), clf, hold on, box on
label_strings = {'(a) 100\% starfish culled at 168 reefs', '(b) 50\% starfish culled at 336 reefs', ...
    '(c) 25\% starfish culled at 672 reefs', '(d) 7.72\% starfish culled at 2175 reefs'};
for i = 1:size(control_plot, 2)
    sp = subplot(2, 2, i); hold on, box on
    % Plot outline of Australia
    pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
    % Plot reef locations by color depending on control effort
    for j = 1:size(control_plot, 1)
        if control_plot(j, i) > 0
            p1 = plot(lon(j), lat(j), '.', 'MarkerSize', 12, 'Color', viridis_palette_4(2, :));
        else
            p2 = plot(lon(j), lat(j), '.', 'MarkerSize', 12, 'Color', viridis_palette_4(3, :));
        end
    end
    % Focus the figure on GBR and QLD
    xlim([140, 155])
    ylim([-26, -8])
    set(gca, 'YTick', -26:6:-8);
    % Add labels
    set(gca, 'BoxStyle', 'Full');
    xlabel(label_strings(i), 'Interpreter', 'Latex', 'Fontsize', axis_FS-2)
    title(['Strategy ', num2str(i)], 'Interpreter', 'Latex', 'Fontsize', title_FS-3)
    % Legend for first plot only
    if i == 1
        [h, icons] = legend([p1 p2], 'Starfish culled', 'No starfish culled', 'Interpreter', 'Latex', 'Fontsize', title_FS-2);
        icons = findobj(icons, 'Type', 'line');
        icons = findobj(icons, 'Marker', 'none', '-xor');
        set(icons, 'MarkerSize', 20)
    end
end
