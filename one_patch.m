clear all

% x = coral cover
% y = COTS

% Define problem parameters
alpha_x = 0.05;     % interaction constant in equation for coral
gamma_y = 0.01;     % ineraction constant in equation for COTS
x_0 = 20;           % initial coral cover
y_0 = 10;           % initial COTS 

% Initialise time vectors
t_0 = 0;
t_start = 1;
t_end = 15;
t_vec = 0:1:15;

% Initialise vecors for x and y
x = zeros(1, t_end+1);
y = zeros(1, t_end+1);
x(1) = x_0;
y(1) = y_0;

% Initialise vector for control effort
k_0 = 0.3;
k = k_0 * ones(1, length(t_vec));

% Loop over and calculate population size through time
for t = t_start:t_end
    % Calculate the population for the next year
    x(t+1) = x(t) - alpha_x * x(t) * y(t);
    y(t+1) = y(t) - gamma_y * x(t) * y(t) - k(t) * y(t);
    
    % Check they aren't less than zero
    if x(t+1) < 0
        x(t+1) = 0;
    end
    if y(t+1) < 0
        y(t+1) = 0;
    end
end

% Plotting
figure(2), clf, hold on
plot(t_vec, x, 'Linewidth', 2)
plot(t_vec, y, 'Linewidth', 2)
title('One-Patch Model')
xlabel('Time (t)')
ylabel('Population Size')
legend('Coral Cover', 'COTS', 'Location', 'NorthEastOutside')