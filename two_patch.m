%% Initialisation
clear all

% Define problem parameters
alpha_x = 0.1;         % interaction constant in equation for coral
alpha_y = 0.01;         % interaction constant in equation for COTS
alpha_S = 1;            % constant in COTS survival term
beta_S = 1;             % constant in COTS survival term      
num_reefs = 2;          % number of reefs      

% Initialise time vector
t_0 = 0;
t_start = 1;
t_end = 20;                                 % time in years
t_vec = t_0:1:t_end;
% Initialise a vector for each reef area
A = [100; 70];                              % reef area

% Initialise matrices for coral and COTS
x = zeros(num_reefs, length(t_vec));        % coral cover
y = zeros(num_reefs, length(t_vec));        % COTS
x_0 = [20; 30];                             % initial coral cover
y_0 = [15; 25];                             % initial COTS 
x(:, 1) = x_0;
y(:, 1) = y_0;

% Initialise matrix for control effort
k_0 = [0.3; 0.3]; 
k = k_0 * ones(1, length(t_vec)-1); % control effort

% Initialise matrices for larval recruitment and survival for coral
r_x = 1;
c_x = [0,   1;                              % coral larval dispersal
       0.5, 0.5];       
R = zeros(num_reefs, length(t_vec));        % coral larval recruitment
phi = zeros(num_reefs, length(t_vec));      % actual larval recruitment

% Initialise matrices for larval recruitment for COTS
r_y = 1;
c_y_0 = [20; 10];
c_y = c_y_0 * ones(1, length(t_vec));       % COTS larval dispersal
S = zeros(num_reefs, length(t_vec));        % COTS larval recruitment


%% Solve
% Loop over and calculate population size through time
for t = t_start:t_end
    % Loop over each reef
    for i = 1:num_reefs
        % Calculate COTS larval recruitment
        sum_S = 0;
        for j = 1:num_reefs
            sum_S = sum_S + c_y(j, i) * y(j, t) * A(j) * r_y;
        end
        S(i, t+1) = 1/A(i) * sum_S;
        
        % Calculate coral larval recuitment
        sum_R = 0;
        for j = 1:num_reefs
            sum_R = sum_R + c_x(j, i) * x(j, t) * A(j) * r_x;
        end
        R(i, t+1) = 1/A(i) * sum_R;        
        
        % Calculate coral settlement based on area left on reef
        if R(i, t+1) <= (A(i) - x(i, t))
            phi(i, t+1) = R(i, t+1);
        else
            phi(i, t+1) = A(i) - x(i, t);
        end
        
        % Decide on control effort 
        if y(i, t) < 5
            k(i, t) = 0;
        end
        
        % Calculate the population for the next year
        x(i, t+1) = x(i, t) * (1 - alpha_x * y(i, t)) + phi(i, t+1);
        y(i, t+1) = y(i, t) * (1 - alpha_y * x(i, t) - k(i, t)) ...
                        + (alpha_S * S(i, t+1))/(1 + beta_S * S(i, t+1));
        
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
    xlabel('Time (in years)', 'Interpreter', 'Latex', 'Fontsize', 12)
    ylabel('Population Size', 'Interpreter', 'Latex', 'Fontsize', 12)
    title(['Reef ', num2str(i)], 'Interpreter', 'Latex', 'Fontsize', 13)
    legend('Coral Cover', 'COTS')
end
