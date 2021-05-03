%% Initialisation
clear all

% Define problem parameters
alpha_x = 1;            % growth rate of coral 
beta_x = 0.01;          % mortality rate of coral from starfish
alpha_y = 0.02;         % starfish response to coral
beta_y = 0.01;          % natural mortality rate of starfish
alpha_S = 1;            % constant in starfish larval survival term
beta_S = 1;             % constant in starfish larval survival term      
num_reefs = 2;          % number of reefs      

% Initialise time vector
t_0 = 0;
t_start = 1;
t_end = 10;                                 % time in years
t_vec = t_0:1:t_end;

% Initialise a vector for each reef area
A = [100; 70];                              % reef area

% Initialise matrices for coral and starfish
x = zeros(num_reefs, length(t_vec));        % coral cover each year
x_end = zeros(num_reefs, length(t_vec)-1);  % coral cover at end of year
y = zeros(num_reefs, length(t_vec));        % starfish each year
y_end = zeros(num_reefs, length(t_vec)-1);  % starfish at end of year
x_0 = [25; 30];                             % initial coral cover
y_0 = [15; 25];                             % initial starfish 
x(:, 1) = x_0;
y(:, 1) = y_0;

% Initialise matrix for control effort - note for now this is percentage
k_0 = [0.2; 0.1]; 
k = k_0 * ones(1, length(t_vec)-1);         % control effort

% Larval recruitment and survival for coral
r_x = 5;                                    % production rate
c_x = [0.2, 0.2;
       0.3, 0.1];                           % coral larval dispersal
% We need to account for the percentage that die by floating off
c_x_dead = 1 - sum(c_x, 2);                     
R = zeros(num_reefs, length(t_vec));        % coral larval recruitment
phi = zeros(num_reefs, length(t_vec));      % actual larval recruitment

% Larval recruitment for starfish
r_y = 3;                                    % production rate
c_y = [0.1, 0.4;
       0.5, 0.1];                           % starfish larval dispersal
% We need to account for the percentage that die by floating off
c_y_dead = 1 - sum(c_y, 2); 
S = zeros(num_reefs, length(t_vec));        % starfish larval recruitment


%% Solve
% Loop over and calculate population size through time
for t = t_start:t_end
    % Loop over each reef and calculate population at the end of last year
    for i = 1:num_reefs
        % Calculate population sizes
        x_end(i, t) = x(i, t) * (alpha_x - beta_x * y(i, t));
        y_end(i, t) = y(i, t) * (alpha_y * x(i, t) - beta_y) - k(i, t) * y(i, t);
        % Round no. of starfish to nearest integer
        y_end(i, t) = round(y_end(i, t));

        % Check coral cover isn't less than zero
        if x_end(i, t) < 0
            x_end(i, t) = 0;
        end
        % Check coral cover isn't larger than the reef size
        if x_end(i, t) > A(i)
            x_end(i, t) = A(i);
        end
        % Check no. of starfish isn't less than zero
        if y_end(i, t) < 0
            y_end(i, t) = 0;
        end
    end
    
    % Loop over each reef and calculate population this year based on last
    % year
    for i = 1:num_reefs
        % Calculate starfish larval recruitment
        for j = 1:num_reefs
            S(i, t+1) = S(i, t+1) + c_y(j, i) * y_end(j, t) * r_y;
        end
       
        % Calculate coral larval recuitment
        for j = 1:num_reefs
            R(i, t+1) = R(i, t+1) + c_x(j, i) * x_end(j, t) * r_x;
        end
        % Calculate coral settlement based on area left on reef
        if R(i, t+1) <= (A(i) - x_end(i, t))
            phi(i, t+1) = R(i, t+1);
        else
            phi(i, t+1) = A(i) - x_end(i, t);
        end
        
        % Calculate the population for the next year
        x(i, t+1) = x_end(i, t) + phi(i, t+1);
        y(i, t+1) = y_end(i, t) + (alpha_S * S(i, t+1))/(1 + beta_S * S(i, t+1));
        % Round number of starfish to nearest integer
        y(i, t+1) = round(y(i, t+1));
        
        % Check coral cover isn't less than zero
        if x(i, t+1) < 0
            x(i, t+1) = 0;
        end
        % Check coral cover isn't larger than the reef size
        if x(i, t+1) > A(i)
            x(i, t+1) = A(i);
        end
        % Check no. of starfish isn't less than zero
        if y(i, t+1) < 0
            y(i, t+1) = 0;
        end
    end
end


%% Plotting
% Predator and prey populations over time
figure(1), clf, hold on
for i = 1:num_reefs
    subplot(1, num_reefs, i), hold on
    plot(t_vec, x(i, :), '-o', 'Linewidth', 1.5)
    plot(t_vec, y(i, :), '-o', 'Linewidth', 1.5)
    xlabel('Time (in years)', 'Interpreter', 'Latex', 'Fontsize', 12)
    ylabel('Population Size', 'Interpreter', 'Latex', 'Fontsize', 12)
    title(['Reef ', num2str(i)], 'Interpreter', 'Latex', 'Fontsize', 13)
    legend('Coral Cover (sqm)', 'Starfish (no.)', 'Location', 'NorthWest')
end

% Predator-prey phase plane
figure(2), clf, hold on
for i = 1:num_reefs
    plot(x(i, :), y(i, :), '-', 'Linewidth', 1)
end
xlabel('Coral cover (sqm)', 'Interpreter', 'Latex', 'Fontsize', 12)
ylabel('Number of starfish', 'Interpreter', 'Latex', 'Fontsize', 12)
title('Predator-Prey Phase Plane', 'Interpreter', 'Latex', 'Fontsize', 13)
legend('Reef 1', 'Reef 2')
xlim([0 max(max(x))+20])
