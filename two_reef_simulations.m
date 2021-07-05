clear all;

%% No Control Effort For 15 Years
% PARAMETERS
num_reefs = 2;                  % number of reefs 
t_end = 15;                     % time in years
control_effort = 0;             % no control effort

% Parameters struct
params.alpha_c = 1.01;          % growth rate of coral 
params.beta_sc = 0.01;          % mortality rate of coral from starfish
params.beta_cs = 0.02;          % starfish response to coral
params.alpha_s = 0.01;          % natural mortality rate of starfish
params.rho = 0.8;               % proliferation rate of starfish larvae
params.K = [10; 12];            % starfish larvae carrying capacity
params.r_c = 1.2;               % coral larvae production rate
params.r_s = 1;                 % starfish larvae production rate
params.kappa_c = [0.2, 0.2;
    0.3, 0.1];                  % coral larvae connectivity matrix
params.kappa_s = [0.1, 0.1;
    0.2, 0.1];                  % starfish larvae connectivity matrix
params.A = [100; 70];           % reef area i.e. coral carrying capacity

% Initial system state
initial_state.x_0 = [25; 30];   % initial coral cover
initial_state.y_0 = [15; 25];   % initial starfish popultion


% SOLVE AND PLOT
[t_vec, x, y] = simulate_reefs(num_reefs, t_end, params, initial_state, control_effort);

% Predator-prey populations over time
figure(1), clf, hold on
for i = 1:num_reefs
    subplot(1, num_reefs, i), hold on
    plot(t_vec, x(i, :), '-o', 'Linewidth', 1.5)
    plot(t_vec, y(i, :), '-o', 'Linewidth', 1.5)
    yline(params.A(i), 'k--')
    xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', 13)
    ylabel('Population size', 'Interpreter', 'Latex', 'Fontsize', 13)
    title(['Reef ', num2str(i), ' without control'], 'Interpreter', ...
        'Latex', 'Fontsize', 14) 
    if i == 1
        ylim([0 160])
    end
    if i == 2
        legend('Coral cover (m^2)', 'No. of starfish', ...
               'Carrying capacity of coral', 'Location', 'NorthEast')
        ylim([0 120]) 
    end
end


%% No Control Effort For 50 Years

% PARAMTERS
% Update any parameters
t_end = 50;

% SOLVE AND PLOT
[t_vec_50, x_50, y_50] = simulate_reefs(num_reefs, t_end, params, initial_state, control_effort);

% Predator-prey phase plane
figure(2), clf, hold on
for i = 1:num_reefs
    subplot(1, num_reefs, i), hold on
    plot(x_50(i, 1), y_50(i, 1), '.', 'Color', '#0072BD', 'MarkerSize', 20)
    plot(x_50(i, :), y_50(i, :), 'Color', '#0072BD', 'Linewidth', 1)
    xlabel('Coral cover ($m^2$)', 'Interpreter', 'Latex', 'Fontsize', 13)
    ylabel('Number of starfish', 'Interpreter', 'Latex', 'Fontsize', 13)
    title(['Reef ', num2str(i)], 'Interpreter', 'Latex', 'Fontsize', 14)
    if i == 2
        legend('Initial Population')
    end
end

%% Constant Control Effort Over 15 Years

% PARAMETERS
% Update any parameters
t_end = 15;
k = [0.2; 0.15];
control_effort = k * ones(1, t_end);

% SOLVE AND PLOT
[t_vec_cc, x_cc, y_cc] = simulate_reefs(num_reefs, t_end, params, initial_state, control_effort);

% PLOT 1
% Predator-prey populations over time for constant control effort
figure(3), clf, hold on
for i = 1:num_reefs
    subplot(1, num_reefs, i), hold on
    plot(t_vec_cc, x_cc(i, :), '-o', 'Linewidth', 1.5)
    plot(t_vec_cc, y_cc(i, :), '-o', 'Linewidth', 1.5)
    yline(params.A(i), 'k--')
    xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', 13)
    ylabel('Population size', 'Interpreter', 'Latex', 'Fontsize', 13)
    title(['Reef ', num2str(i), ' with ' num2str(k(i)*100), '\% control'], ... 
        'Interpreter', 'Latex', 'Fontsize', 14) 
    if i == 1
        ylim([0 140])
    end
    if i == 2
        legend('Coral cover (m^2)', 'No. of starfish', ...
               'Carrying capacity of coral', 'Location', 'NorthEast')
        ylim([0 120]) 
    end
end


% PLOT 2
% Control effort comparison for reef 1 
figure(4), clf, hold on
i = 1;
% No control effort
subplot(1, 2, 1), hold on
plot(t_vec, x(i, :), '-o', 'Linewidth', 1.5)
plot(t_vec, y(i, :), '-o', 'Linewidth', 1.5)
yline(params.A(i), 'k--')
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', 13)
ylabel('Population size', 'Interpreter', 'Latex', 'Fontsize', 13)
title(['Reef ', num2str(i), ' without control'], 'Interpreter', ...
        'Latex', 'Fontsize', 14) 
ylim([0 160])
% Constant control effort
subplot(1, 2, 2), hold on
plot(t_vec_cc, x_cc(i, :), '-o', 'Linewidth', 1.5)
plot(t_vec_cc, y_cc(i, :), '-o', 'Linewidth', 1.5)
yline(params.A(i), 'k--')
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', 13)
ylabel('Population size', 'Interpreter', 'Latex', 'Fontsize', 13)
title(['Reef ', num2str(i), ' with ' num2str(k(i)*100), '\% control'], ...
        'Interpreter', 'Latex', 'Fontsize', 14)     
legend('Coral cover (m^2)', 'No. of starfish', ...
    'Carrying capacity of coral', 'Location', 'NorthEast')
ylim([0 140])


% PLOT 3
% Control effort comparison for reef 2
figure(5), clf, hold on
i = 2;
% No control effort
subplot(1, 2, 1), hold on
plot(t_vec, x(i, :), '-o', 'Linewidth', 1.5)
plot(t_vec, y(i, :), '-o', 'Linewidth', 1.5)
yline(params.A(i), 'k--')
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', 13)
ylabel('Population size', 'Interpreter', 'Latex', 'Fontsize', 13)
title(['Reef ', num2str(i), ' without control'], 'Interpreter', ...
        'Latex', 'Fontsize', 14) 
ylim([0 120])
% Constant control effort
subplot(1, 2, 2), hold on
plot(t_vec_cc, x_cc(i, :), '-o', 'Linewidth', 1.5)
plot(t_vec_cc, y_cc(i, :), '-o', 'Linewidth', 1.5)
yline(params.A(i), 'k--')
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', 13)
ylabel('Population size', 'Interpreter', 'Latex', 'Fontsize', 13)
title(['Reef ', num2str(i), ' with ' num2str(k(i)*100), '\% control'], ...
        'Interpreter', 'Latex', 'Fontsize', 14)  
legend('Coral cover (m^2)', 'No. of starfish', ...
    'Carrying capacity of coral', 'Location', 'NorthEast')
ylim([0 120])