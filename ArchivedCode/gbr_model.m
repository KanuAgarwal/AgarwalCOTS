clear all

%% PARAMETERISATION

% Get connectivity matrices
load IdentifyKeySources/ConnectivityMatrices_Model_A_2002_P7

% PARAMETERS
num_reefs = 2175;               % number of reefs 
t_end = 25;                     % time in years
control_effort = 0;             % no control effort

% Parameters struct: store all parameters in a struct
% Constant parameter values for all reefs
params.alpha_c = 1.01;          % growth rate of coral 
params.beta_sc = 0.01;          % mortality rate of coral from starfish
params.beta_cs = 0.02;          % starfish response to coral
params.alpha_s = 0.01;          % natural mortality rate of starfish
params.rho = 0.8;               % proliferation rate of starfish larvae
params.r_c = 1.2;               % coral larvae production rate
params.r_s = 1;                 % starfish larvae production rate

% Reef dependent variables
% starfish larvae carrying capacity
params.K = 15 * ones(num_reefs, 1);         % Need to find
% coral larvae connectivity matrix
params.kappa_c = psurv_d02_1122_P7;         % From Bode et al. (2012)
% starfish larvae connectivity matrix
params.kappa_s = psurv_d02_1122_P7;         % From Bode et al. (2012)
% reef area i.e. coral carrying capacity
params.A = 100 * ones(num_reefs, 1);        % Need to find

% Initial system state
% initial coral cover - need to input from data
initial_state.x_0 = 30 * ones(num_reefs, 1);
% initial starfish popultion - need to input from data
initial_state.y_0 = 20 * ones(num_reefs, 1);   


%% SOLVE
[t_vec, x, y] = simulate_reefs(num_reefs, t_end, params, initial_state, control_effort);


%% PLOTS

% Figure 1: Coral heatmap
figure(1), clf, hold on
imagesc(x)
title('Coral populations')
xlabel('Time (years)')
ylabel('Reef')
colorbar
xlim([0 26])
ylim([0 num_reefs])

% Figure 2: Starfish heatmap
figure(2), clf, hold on
imagesc(y)
title('Starfish populations')
xlabel('Time (years)')
ylabel('Reef')
colorbar
xlim([0 26])
ylim([0 num_reefs])


% Figure 3: GBR animation with coral populations over time
figure(3), clf, hold on

% Load in australian map outline
load AustOutline

% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);

% Load in lat long data
load IdentifyKeySources/original_centroids

% Plot reef locations in grey
plot(lg, lt, '.', 'Markersize', 10, 'Color', 0.7.*ones(1,3))

% % Load in the reef outline data
% load ReefOutline
% 
% % Plot the reefs on the GBR as blue outlines
% plot(ReefRaw(:, 1), ReefRaw(:, 2), 'b')

% Focus the figure on GBR and QLD
xlim([140, 155])
ylim([-26, -8])
