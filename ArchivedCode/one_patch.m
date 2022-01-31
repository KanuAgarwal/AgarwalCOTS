clear all

% x = coral cover
% y = COTS

% Define problem parameters
alpha_x = 0.05;     % interaction constant in equation for coral
alpha_y = 0.01;     % interaction constant in equation for COTS
% alpha = 1;
% beta = 1;
% delta = 1;
% gamma = 1;
x_0 = 20;           % initial coral cover
y_0 = 10;           % initial COTS 

% Initialise time vector
t_0 = 0;
t_start = 1;
t_end = 15;
t_vec = t_0:1:t_end;

% Initialise vectors for x and y
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
    x(t+1) = x(t) * (1 - alpha_x * y(t));
    y(t+1) = y(t) * (1 - alpha_y * x(t) - k(t));    
%     x(t+1) = x(t) * (alpha - beta * y(t));
%     y(t+1) = y(t) * (delta * x(t) - gamma - k(t));
    
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
plot(t_vec, x, '-o', 'Linewidth', 1.5)
plot(t_vec, y, '-o', 'Linewidth', 1.5)
title('One-Patch Model')
xlabel('Time (t)')
ylabel('Population Size')
legend('Coral Cover', 'COTS', 'Location', 'NorthEastOutside')