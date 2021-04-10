%% Initialisation
clear all

% Define problem parameters
alpha_x = 0.05;         % interaction constant in equation for coral
alpha_y = 0.01;         % interaction constant in equation for COTS
alpha_s = 1;            % constant in COTS survival term
beta_s = 1;             % constant in COTS survival term      
num_reefs = 2;          % number of reefs      

% Initialise time vector
t_0 = 0;
t_start = 1;
t_end = 20;
t_vec = t_0:1:t_end;

% Initialise a vector for each reef area
A_0 = [100; 70];
A = A_0 * ones(1, length(t_vec));              % reef area

% Initialise matrices for coral and COTS
x = zeros(num_reefs, length(t_vec));        % coral cover
y = zeros(num_reefs, length(t_vec));        % COTS
x_0 = [20; 30];                             % initial coral cover
y_0 = [15; 25];                             % initial COTS 
x(:, 1) = x_0;
y(:, 1) = y_0;

% Initialise matrix for control effort
k_0 = 0.5; 
k = k_0 * ones(num_reefs, length(t_vec)-1);   % control effort

% Initialise matrices for larval recruitment and survival
R_0 = [3; 2];
R = R_0 * ones(1, length(t_vec));           % coral larval recruitment
S_0 = [10; 4];
S = S_0 * ones(1, length(t_vec));           % COTS larval recruitment
phi = zeros(num_reefs, length(t_vec));      % coral larval survival 


%% Solve
% Loop over and calculate population size through time
for t = t_start:t_end
    % Loop over each reef
    for i = 1:num_reefs
        % Calculate coral settlement based on area left on reef
        if R(i, t) <= (A(i) - x(i, t))
            phi(i, t+1) = R(i, t);
        else
            phi(i, t+1) = A(i) - x(i, t);
        end
        
        % Decide on control effort 
        if y(i, t) < 10
            k(i, t) = 0.3;
        elseif y(i, t) < 5
            k(i, t) = 0;
        end
        
        % Calculate the population for the next year
        x(i, t+1) = x(i, t) * (1 - alpha_x * y(i, t)) + phi(i, t+1);
        y(i, t+1) = y(i, t) * (1 - alpha_y * x(i, t) - k(i, t)) ...
                        + (alpha_s * S(i, t+1))/(1 + beta_s * S(i, t+1));
        
        % Check they aren't less than zero
        if x(i, t+1) < 0
            x(i, t+1) = 0;
        end
        if y(i, t+1) < 0
            y(i, t+1) = 0;
        end
    end
end


%% Plotting
figure(3), clf, hold on
for i = 1:num_reefs
    subplot(1, num_reefs, i), hold on
    plot(t_vec, x(i, :), '-o', 'Linewidth', 1.5)
    plot(t_vec, y(i, :), '-o', 'Linewidth', 1.5)
    ylim([0 40])
    xlabel('Time (in years)', 'Interpreter', 'Latex', 'Fontsize', 12)
    ylabel('Population Size', 'Interpreter', 'Latex', 'Fontsize', 12)
    title(['Reef ', num2str(i)], 'Interpreter', 'Latex', 'Fontsize', 13)
    legend('Coral Cover', 'COTS')
end
