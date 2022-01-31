%% Parameter values and initial values
clear all

% Define problem parameters
alpha_x = 1;            % growth rate of coral 
beta_x = 0.01;          % mortality rate of coral from starfish
alpha_y = 0.02;         % starfish response to coral
beta_y = 1;             % natural mortality rate of starfish
alpha_S = 1;            % constant in starfish larval survival term
beta_S = 1;             % constant in starfish larval survival term  

% Turn into params struct
params.alpha_x = alpha_x;
params.beta_x = beta_x;
params.alpha_y = alpha_y;
params.beta_y = beta_y;
params.alpha_S = alpha_S;
params.beta_S = beta_S;

% Initialise time vector
t_0 = 0;
t_end = 10;                                 % time in years
t_span = [t_0 t_end];                       % Span of time

% Number of reefs in problem
num_reefs = 2;                              % number of reefs      

% Initialise a vector for each reef area
A = [100; 70];                              % reef area

% Initial coral cover and starfish at reef 1
xy_1_0 = [25; 15];
% Initial coral cover and starfish at reef 2
xy_2_0 = [30; 25];


%% Control effort and larval dispersal
% Initialise matrix for control effort - note for now this is percentage
k_0 = [0; 0]; 

% Larval recruitment and survival for coral
r_x = 5;                                    % production rate
c_x = [0.2, 0.2;
       0.3, 0.1];                           % coral larval dispersal
% We need to account for the percentage that die by floating off
c_x_dead = 1 - sum(c_x, 2);                     

% Larval recruitment for starfish
r_y = 3;                                    % production rate
c_y = [0.1, 0.4;
       0.5, 0.1];                           % starfish larval dispersal
% We need to account for the percentage that die by floating off
c_y_dead = 1 - sum(c_y, 2); 


%% Solving using ode45
% Reef 1
[t_vals_1, xy_1] = ode45(@(t, xy) ode_model(t, xy, params, 0.5, 0, 0), ...
                       t_span, xy_1_0);
x_vals_1 = xy_1(:, 1);
y_vals_1 = xy_1(:, 2);

% Reef 2
[t_vals_2, xy_2] = ode45(@(t, xy) ode_model(t, xy, params, 0.2, 0, 0), ...
                       t_span, xy_2_0);
x_vals_2 = xy_2(:, 1);
y_vals_2 = xy_2(:, 2);


%% Plotting
% Predator and prey populations over time
figure(1), clf, hold on
% Reef 1
subplot(1, num_reefs, 1), hold on
plot(t_vals_1, x_vals_1, 'Linewidth', 2)
plot(t_vals_1, y_vals_1, 'Linewidth', 2)
yline(A(1), 'k--')
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', 12)
ylabel('Population size', 'Interpreter', 'Latex', 'Fontsize', 12)
title(['Reef ', num2str(1)], 'Interpreter', 'Latex', 'Fontsize', 13) 
legend('Coral cover (m^2)', 'No. of starfish', 'Location', 'NorthEast')
% Reef 2
subplot(1, num_reefs, 2), hold on
plot(t_vals_2, x_vals_2, 'Linewidth', 2)
plot(t_vals_2, y_vals_2, 'Linewidth', 2)
yline(A(2), 'k--')
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', 12)
ylabel('Population size', 'Interpreter', 'Latex', 'Fontsize', 12)
title(['Reef ', num2str(2)], 'Interpreter', 'Latex', 'Fontsize', 13)

% Predator-prey phase plane
figure(2), clf, hold on
% Reef 1
subplot(1, num_reefs, 1), hold on
plot(x_vals_1, y_vals_1, 'Linewidth', 2)
xlabel('Coral cover ($m^2$)', 'Interpreter', 'Latex', 'Fontsize', 12)
ylabel('Number of starfish', 'Interpreter', 'Latex', 'Fontsize', 12)
title(['Reef ', num2str(1)], 'Interpreter', 'Latex', 'Fontsize', 13)
% Reef 2
subplot(1, num_reefs, 2), hold on
plot(x_vals_2, y_vals_2, 'Linewidth', 2)
xlabel('Coral cover ($m^2$)', 'Interpreter', 'Latex', 'Fontsize', 12)
ylabel('Number of starfish', 'Interpreter', 'Latex', 'Fontsize', 12)
title(['Reef ', num2str(2)], 'Interpreter', 'Latex', 'Fontsize', 13)


%% ODE functions

function xy_dash = ode_model(t, xy, params, k, phi, S)
% Parameters
alpha_x = params.alpha_x;
beta_x = params.beta_x;
alpha_y = params.alpha_y;
beta_y = params.beta_y;
alpha_S = params.alpha_S;
beta_S = params.beta_S;

% x and y
x = xy(1);
y = xy(2);

% ODE 
xy_dash = [x * (alpha_x - beta_x * y) + phi; 
           y * (alpha_y * x - beta_y - k) + (alpha_S * S)/(1 + beta_S * S)];
end