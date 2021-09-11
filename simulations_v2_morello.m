clear all

% PARAMETERISATION ========================================================
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
t_end = 25;                     % time in years

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
params.p_2_f = 10;              % effect of COTS on fast-growing coral
% params.r_m = 0.1;               % intrinsic growth rate of slow-growing coral
% params.K_m = 500;               % carrying capacity of slow-growing coral
% params.p_2_m = 8;               % effect of COTS on slow-growing coral

% Known or arbitrarily chosen by me
params.r_c = 0.1;               % coral larvae reproduction rate
params.r_s = 0;                 % starfish larvae reproduction rate
% params.rho = 1;                 % proliferation rate of starfish larvae 
% params.K = 1;                   % starfish larvae carrying capacity

% Connectivity matrices from Bode et al. (2012)
params.omega_c = psurv_d02_1122_P7;             % coral
params.omega_s = psurv_d02_1122_P7;             % starfish

% Individual reef parameters 
% Get number of reefs from connectivity matrices
num_reefs = length(params.omega_c);
% % Area of each reef i.e. coral carrying capacity
% params.A = 1;


% INITIAL SYSTEM STATE ----------------------------------------------------
% CORAL
% Percentage of fast-growing coral = carrying capacity
initial_state.C_0_f = 0.8 * params.K_f * ones(num_reefs, 1);
% % Biomass of slow-growing coral = carrying capacity
% initial_state.C_0_m = params.K_m * ones(num_reefs, 1);  

% STARFISH
% Number of COTS aged 2+ = estimated
% initial_state.N_0_2 = 0.505 * ones(num_reefs, 1);
initial_state.N_0_2 = zeros(num_reefs, 1);

% % Start the outbreak
% % Find index of reef with the most connectivity
% [~, index] = max(sum(params.omega_c, 2));
% % Put some starfish there
% initial_state.N_0_2(index, 1) = 200;
% % Pick other reefd to put starfish at
% index_2 = 304;
% initial_state.N_0_2(index_2, 1) = 200;
% index_3 = 1004;
% initial_state.N_0_2(index_3, 1) = 200;

% Number of COTS aged 1 = no. of age 2+ COTS * function
% initial_state.N_0_1 = initial_state.N_0_2 * exp(params.M_cots);
initial_state.N_0_1 = zeros(num_reefs, 1);
% Number of COTS aged 0 = no. of age 2+ COTS * function
% initial_state.N_0_0 = initial_state.N_0_2 * exp(2*params.M_cots);
initial_state.N_0_0 = zeros(num_reefs, 1);


% CONTROL EFFORT ----------------------------------------------------------
control_effort = 0;             % no control effort



% SOLVE ===================================================================
% Solve using function which runs simulations
[t_vec, C_y_f, N_y_2, N_y_1, N_y_0] = simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort);

% % Calculate coral cover from coral biomass
% C_y_f_cover = C_y_f.^(2/3);
% C_y_m_cover = C_y_m.^(2/3);

% Calculate coral cover and cots over time
coral_over_time = sum(C_y_f, 1);
starfish_over_time = sum(N_y_2, 1);


% PLOTS ===================================================================
% % Coral biomass -----------------------------------------------------------
% figure, clf, hold on, grid on
% plot(t_vec, C_y_f(1, :), 'Linewidth', 2)
% plot(t_vec, C_y_m(1, :), 'Linewidth', 2)
% xlabel('Time (years)')
% ylabel('Coral biomass')
% title('Coral biomass at reef 1')
% legend('Fast-growing coral', 'Slow-growing coral')

% % Fast-growing coral cover ------------------------------------------------
% figure, clf, hold on, grid on
% plot(t_vec, C_y_f_cover(1, :), 'Linewidth', 2)
% xlabel('Time (years)')
% ylabel('Coral cover')
% title('Fast-growing coral cover at reef 1')

% % Slow-growing coral cover ------------------------------------------------
% figure, clf, hold on, grid on
% plot(t_vec, C_y_m_cover(1, :), 'Linewidth', 2, 'Color', [0.8500 0.3250 0.0980])
% xlabel('Time (years)')
% ylabel('Coral cover')
% title('Slow-growing coral cover at reef 1')

% % All COTS ----------------------------------------------------------------
% figure, clf, hold on, grid on
% plot(t_vec, N_y_2(1, :), 'Linewidth', 2)
% plot(t_vec, N_y_1(1, :), 'Linewidth', 2)
% plot(t_vec, N_y_0(1, :), 'Linewidth', 2)
% xlabel('Time (years)')
% ylabel('Number of COTS')
% legend('Age 2+ COTS', 'Age 1 COTS', 'Age 0 COTS')
% title('COTS population at reef 1')

% % Age 2 COTS --------------------------------------------------------------
% figure, clf, hold on, grid on
% plot(t_vec, N_y_2(1, :), 'Linewidth', 2)
% xlabel('Time (years)')
% ylabel('Number of COTS')
% title('Age 2+ COTS population at reef 1')

% Fast-growing coral heatmap ----------------------------------------------
figure(1), clf, hold on
imagesc(C_y_f)
title('Fast-growing coral biomass')
xlabel('Time (years)')
ylabel('Reef')
colorbar
xlim([1 t_end+1])
ylim([0 num_reefs])

% Age 2+ COTS heatmap -----------------------------------------------------
figure(2), clf, hold on
imagesc(N_y_2)
title('Age 2+ COTS')
xlabel('Time (years)')
ylabel('Reef')
colorbar
xlim([1 t_end+1])
ylim([0 num_reefs])

% Total coral cover over time ---------------------------------------------
figure(3), clf, hold on, grid on
plot(t_vec, coral_over_time, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total coral cover')
title('Total coral cover on GBR over time')

% Total starfish over time ------------------------------------------------
figure(4), clf, hold on, grid on
plot(t_vec, starfish_over_time, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total starfish')
title('Total Age 2+ starfish on GBR over time')

% % Outbreak initiation on GBR ----------------------------------------------
% figure(5), clf, hold on, box on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
% % Plot reef locations in grey
% plot(lat, lon, '.', 'Markersize', 10, 'Color', 0.7.*ones(1,3))
% % Reef where outbreak is starting in this simulation
% plot(lat(index), lon(index), 'b.', 'Markersize', 10)
% plot(lat(index_2), lon(index_2), 'b.', 'Markersize', 10)
% plot(lat(index_3), lon(index_3), 'b.', 'Markersize', 10)
% % Focus the figure on GBR and QLD
% xlim([140, 155])
% ylim([-26, -8])
% title('Starfish outbreak locations')

% % Coral cover on GBR ------------------------------------------------------
% figure(6), clf, hold on, box on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
% title('Coral cover on GBR after 15 years')
% xlim([140, 155])
% ylim([-26, -8])
% % Plot reef locations by colour based on coral presence
% % for t = 1:length(t_vec)
% for i = 1:length(lat)
%     if C_y_f(i, end) > 0.9
%         plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.4660 0.6740 0.1880])
%     elseif C_y_f(i, end) < 0.01
%         plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.8500 0.3250 0.0980])
%     else
%         plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.9290 0.6940 0.1250])
%     end
% end
% % end

% % Starfish population on GBR ----------------------------------------------
% figure(7), clf, hold on, box on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
% title('Starfish population on GBR after 15 years')
% xlim([140, 155])
% ylim([-26, -8])
% % Plot reef locations by colour based on starfish presence
% % for t = 1:length(t_vec)
% for i = 1:length(lat)
%     if N_y_2(i, end) > 0.5
%         plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.8500 0.3250 0.0980])
%     else
%         plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.4660 0.6740 0.1880])
%     end
% end
% % end


% % Connectivity of reefs ---------------------------------------------------
% figure(8), clf, hold on, grid on
% plot(1:1:2175, sum(params.omega_c, 2))
