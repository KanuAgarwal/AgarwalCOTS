clear all

% PARAMETERISATION
% DATA
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

% PARAMETERS
t_end = 15;                     % time in years
control_effort = 0;             % no control effort

% Parameters struct: store all parameters in a struct
% Constant parameter values for all reefs
% Estimated by Morello et al. (2014)
params.p_tilde = 0.258;         % effect of fast-growing coral on COTS     
params.M_cots = 2.56;           % natural mortality of COTS 
params.p_1_f = 0.129;           % effect of COTS on fast-growing coral
params.p_1_m = 0.268;           % effect of COTS on slow-growing coral

% Known or arbitrarily chosen by Morello et al. (2014)
params.r_f = 0.5;               % intrinsic growth rate of fast-growing coral
params.r_m = 0.1;               % intrinsic growth rate of slow-growing coral
params.K_f = 1;              % carrying capacity of fast-growing coral
params.K_m = 500;               % carrying capacity of slow-growing coral
params.p_2_f = 10;              % effect of COTS on fast-growing coral
params.p_2_m = 8;               % effect of COTS on slow-growing coral

% Known or arbitrarily chosen by me
% params.rho = 1;                 % proliferation rate of starfish larvae 
% params.K = 1;                   % starfish larvae carrying capacity

% Estimated or taken from other papers
% Coral larvae production rate - from Practchett et al. (2019)
params.r_c = 0.1;  
% Starfish larvae production rate - from Pratchett et al. (2021)
params.r_s = 6730;                 

% Reef dependent variables
% Coral larvae connectivity matrix - from Bode et al. (2012)
params.omega_c = psurv_d02_1122_P7;
% Starfish larvae connectivity matrix - from Bode et al. (2013)
params.omega_s = psurv_d02_1122_P7;
% Number of reefs - extract from connectivity matrices
num_reefs = length(params.omega_c);
% params.A = 1;                   % reef area i.e. coral carrying capacity


% INITIAL SYSTEM STATE
% Coral
% biomass of fast-growing coral = carrying capacity
initial_state.C_0_f = 0.8* params.K_f * ones(num_reefs, 1);
% biomass of slow-growing coral = carrying capacity
initial_state.C_0_m = params.K_m * ones(num_reefs, 1);  

% Starfish
% number of COTS aged 2+ = estimated
% initial_state.N_0_2 = 0.505 * ones(num_reefs, 1);
initial_state.N_0_2 = zeros(num_reefs, 1);
% Find index of reef at top
% [~, index] = max(lat);
% % Put some starfish there
% initial_state.N_0_2(492, 1) = 1000;

% number of COTS aged 1 = no. of age 2+ COTS * function
% initial_state.N_0_1 = initial_state.N_0_2 * exp(params.M_cots);
initial_state.N_0_1 = zeros(num_reefs, 1);
% number of COTS aged 0 = no. of age 2+ COTS * function
% initial_state.N_0_0 = initial_state.N_0_2 * exp(2*params.M_cots);
initial_state.N_0_0 = zeros(num_reefs, 1);


% SOLVE
[t_vec, C_y_f, N_y_2, N_y_1, N_y_0] = simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort);

% Calculate coral cover from coral biomass
% C_y_f_cover = C_y_f.^(2/3);
% C_y_m_cover = C_y_m.^(2/3);

% Calculate coral cover and cots over time
coral_over_time = sum(C_y_f, 1);
starfish_over_time = sum(N_y_2, 1);

% PLOTS
% Coral biomass
% figure, clf, hold on, grid on
% plot(t_vec, C_y_f(1, :), 'Linewidth', 2)
% plot(t_vec, C_y_m(1, :), 'Linewidth', 2)
% xlabel('Time (years)')
% ylabel('Coral biomass')
% title('Coral biomass at reef 1')
% legend('Fast-growing coral', 'Slow-growing coral')

% Fast-growing coral cover
% figure, clf, hold on, grid on
% plot(t_vec, C_y_f_cover(1, :), 'Linewidth', 2)
% xlabel('Time (years)')
% ylabel('Coral cover')
% title('Fast-growing coral cover at reef 1')

% Slow-growing coral cover
% figure, clf, hold on, grid on
% plot(t_vec, C_y_m_cover(1, :), 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980])
% xlabel('Time (years)')
% ylabel('Coral cover')
% title('Slow-growing coral cover at reef 1')

% All COTS
% figure, clf, hold on, grid on
% plot(t_vec, N_y_2(1, :), 'Linewidth', 2)
% plot(t_vec, N_y_1(1, :), 'Linewidth', 2)
% plot(t_vec, N_y_0(1, :), 'Linewidth', 2)
% xlabel('Time (years)')
% ylabel('Number of COTS')
% legend('Age 2+ COTS', 'Age 1 COTS', 'Age 0 COTS')
% title('COTS population at reef 1')

% Age 2 COTS 
% figure, clf, hold on, grid on
% plot(t_vec, N_y_2(1, :), 'Linewidth', 2)
% xlabel('Time (years)')
% ylabel('Number of COTS')
% title('Age 2+ COTS population at reef 1')

% Fast-growing coral heatmap
figure(1), clf, hold on
imagesc(C_y_f)
title('Fast-growing coral biomass')
xlabel('Time (years)')
ylabel('Reef')
colorbar
xlim([1 16])
ylim([0 num_reefs])

% Age 2+ COTS heatmap
figure(2), clf, hold on
imagesc(N_y_2)
title('Age 2+ COTS')
xlabel('Time (years)')
ylabel('Reef')
colorbar
xlim([1 16])
ylim([0 num_reefs])

% Total coral cover over time
figure(3), clf, hold on, grid on
plot(t_vec, coral_over_time, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total coral cover')
title('Total coral cover on GBR over time')

% GBR
figure(4), clf, hold on, box on
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
% Plot reef locations in grey
plot(lat, lon, '.', 'Markersize', 10, 'Color', 0.7.*ones(1,3))
% Focus the figure on GBR and QLD
xlim([140, 155])
ylim([-26, -8])
