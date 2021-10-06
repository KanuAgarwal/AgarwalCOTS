% INITIAL SIMULATIONS: Reproduce results from Morello et al. (2014) 

clear all

% PARAMETERISATION ========================================================
% PARAMETERS --------------------------------------------------------------
% How long do we want to run the simulation for
t_end = 17;                     % time in years

% This is for a single reef population
num_reefs = 1;

% Parameters struct: store all parameters in a struct
% Constant parameter values for all reefs
% Estimated by Morello et al. (2014)
params.p_tilde = 0.258;         % effect of fast-growing coral on COTS     
params.M_cots = 2.56;           % natural mortality of COTS 
params.p_1_f = 0.129;           % effect of COTS on fast-growing coral
% params.p_1_m = 0.268;           % effect of COTS on slow-growing coral

% Known or arbitrarily chosen by Morello et al. (2014)
params.r_f = 0.5;               % intrinsic growth rate of fast-growing coral
params.K_f = 1;                 % carrying capacity of fast-growing coral
% params.K_f = 2500;              % carrying capacity of fast-growing coral
params.p_2_f = 10;              % effect of COTS on fast-growing coral
% params.r_m = 0.1;               % intrinsic growth rate of slow-growing coral
% params.K_m = 500;               % carrying capacity of slow-growing coral
% params.p_2_m = 8;               % effect of COTS on slow-growing coral

% Below values not required - but need to enter for function to run
params.r_c = 0;                 % coral larvae reproduction rate
params.r_s = 0;                 % starfish larvae reproduction rate
params.omega_c = 0;             % coral connectivity matric
params.omega_s = 0;             % starfish connectivity matrix


% INITIAL SYSTEM STATE ----------------------------------------------------
% CORAL
% Percentage of fast-growing coral = carrying capacity
initial_state.C_0_f = params.K_f;

% STARFISH
% Number of COTS aged 2+ = estimated
initial_state.N_0_2 = 0.505;

% Initialise age 1 and age 0 COTS based on Morello initial conditions
initial_state.N_0_1 = initial_state.N_0_2 * exp(params.M_cots);
initial_state.N_0_0 = initial_state.N_0_2 * exp(2*params.M_cots);


% CONTROL EFFORT ----------------------------------------------------------
% No control effort
control_effort = 0;


% SOLVE ===================================================================
% Solve using function which runs simulations
[t_vec, C_y_f, N_y_2, N_y_1, N_y_0] = simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort);

% % Calculate coral cover from coral biomass
% C_y_f_cover = C_y_f.^(2/3);


% PLOTS ===================================================================
% Setup -------------------------------------------------------------------
% Create new time vector for ploting
t_plot = 1994:1:2011;

% Fontsizes for plotting
axis_FS = 15;
title_FS = 17;
legend_FS = 12;
ticks_FS = 12;


% % Coral biomass -----------------------------------------------------------
% figure(1), clf, hold on, grid on
% plot(t_plot, C_y_f(1, :), 'Linewidth', 2)
% plot(t_plot, C_y_m(1, :), 'Linewidth', 2)
% xlabel('Time (years)')
% ylabel('Coral biomass')
% title('Coral biomass at reef 1')
% legend('Fast-growing coral', 'Slow-growing coral')

% Fast-growing coral cover ------------------------------------------------
figure(2), clf, hold on, grid on
plot(t_plot, C_y_f(1, :), 'Linewidth', 2)
xlim([1994 2011])
ylim([0 1])
set(gca, 'FontSize', ticks_FS);
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Coral cover (\% of reef area)', 'Interpreter', 'Latex', ...
    'Fontsize', axis_FS)
title('Coral cover at single reef', 'Interpreter', 'Latex', 'Fontsize', title_FS)

% % Slow-growing coral cover ------------------------------------------------
% figure(3), clf, hold on, grid on
% plot(t_plot, C_y_m_cover(1, :), 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980])
% xlabel('Time (years)')
% ylabel('Coral cover')
% title('Slow-growing coral cover at reef 1')

% All COTS ----------------------------------------------------------------
figure(4), clf, hold on, grid on
plot(t_plot, N_y_2(1, :), 'Linewidth', 2)
plot(t_plot, N_y_1(1, :), 'Linewidth', 2)
plot(t_plot, N_y_0(1, :), 'Linewidth', 2)
xlim([1994 2011])
set(gca, 'FontSize', ticks_FS);
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('No. of starfish', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
title('Starfish population at single reef', 'Interpreter', 'Latex', ...
    'Fontsize', title_FS)
legend('Age 2+ (adult)', 'Age 1 (juvenile)', 'Age 0 (larvae)', ...
    'Interpreter', 'Latex', 'Fontsize', legend_FS, 'Location', 'NorthEastOutside')

% Age 2 COTS --------------------------------------------------------------
figure(5), clf, hold on, grid on
plot(t_plot, N_y_2(1, :), 'Linewidth', 2)
xlim([1994 2011])
set(gca, 'FontSize', ticks_FS);
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('No. of age 2+ (adult) starfish', 'Interpreter', 'Latex', ...
    'Fontsize', axis_FS)
title('Adult starfish population at single reef', 'Interpreter', 'Latex', ...
    'Fontsize', title_FS)

