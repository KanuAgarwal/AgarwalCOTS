%% Initialisation
clear all

% Define problem parameters
alpha_c = 1.01;         % growth rate of coral 
beta_sc = 0.01;         % mortality rate of coral from starfish
beta_cs = 0.02;         % starfish response to coral
alpha_s = 0.01;         % natural mortality rate of starfish
rho = 3;                % proliferation rate of starfish larvae 
num_reefs = 2;          % number of reefs      

% Initialise time vector
t_0 = 0;
t_start = 1;
t_end = 15;             % time in years
t_vec = t_0:1:t_end;

% Initialise vectors for carrying capacities
A = [100; 70];          % reef area
K = [2; 2];             % starfish larvae

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
k_0 = [0.4; 0.2]; 
k = k_0 * ones(1, length(t_vec)-1);         % control effort

% Larval recruitment and survival for coral
r_c = 1.2;                                    % production rate
kappa_c = [0.2, 0.2;
           0.3, 0.1];                       % coral larval dispersal
% We need to account for the percentage that die by floating off
kappa_c_dead = 1 - sum(kappa_c, 2);                     
sigma = zeros(num_reefs, length(t_vec));    % coral larval recruitment
phi = zeros(num_reefs, length(t_vec));      % actual larval recruitment

% Larval recruitment for starfish
r_s = 1.5;                                    % production rate
kappa_s = [0.1, 0.1;
           0.2, 0.1];                       % starfish larval dispersal
% We need to account for the percentage that die by floating off
kappa_s_dead = 1 - sum(kappa_s, 2); 
tau = zeros(num_reefs, length(t_vec));      % starfish larval recruitment
gamma = zeros(num_reefs, length(t_vec));    % actual larval recruitment


%% Solve
% Loop over and calculate population size through time
for t = t_start:t_end
    % Loop over each reef and calculate population at the end of last year
    for i = 1:num_reefs
        % Calculate population sizes
        x_end(i, t) = x(i, t) * (alpha_c - beta_sc * y(i, t));
        y_end(i, t) = y(i, t) * (beta_cs * x(i, t) - alpha_s) - k(i, t) * y(i, t);
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
            tau(i, t+1) = tau(i, t+1) + kappa_s(j, i) * y_end(j, t) * r_s;
        end
        % Use Beverton-Holt model for starfish survival
        gamma(i, t+1) = (rho * tau(i, t+1)) / (1 + (rho - 1)/K(i) * tau(i, t+1));
       
        % Calculate coral larval recuitment
        for j = 1:num_reefs
            sigma(i, t+1) = sigma(i, t+1) + kappa_c(j, i) * x_end(j, t) * r_c;
        end
        % Calculate coral settlement based on area left on reef
        if sigma(i, t+1) <= (A(i) - x_end(i, t))
            phi(i, t+1) = sigma(i, t+1);
        else
            phi(i, t+1) = A(i) - x_end(i, t);
        end
        
        % Calculate the population for the next year
        x(i, t+1) = x_end(i, t) + phi(i, t+1);
        y(i, t+1) = y_end(i, t) + gamma(i, t+1);
        
        % Round number of starfish to nearest integer
        if y(i, t+1) >= 1
            y(i, t+1) = round(y(i, t+1));
        elseif y(i, t+1) > 0
            y(i, t+1) = ceil(y(i, t+1));
        % Check it's not less than zero
        elseif y(i, t+1) < 0
            y(i, t+1) = 0;
        end 
        
        % Check coral cover isn't less than zero
        if x(i, t+1) < 0
            x(i, t+1) = 0;
        end
        % Check coral cover isn't larger than the reef size
        if x(i, t+1) > A(i)
            x(i, t+1) = A(i);
        end
    end
end


%% Plotting - No Control Effort
% Predator and prey populations over time
% figure(1), clf, hold on
% for i = 1:num_reefs
%     subplot(1, num_reefs, i), hold on
%     plot(t_vec, x(i, :), '-o', 'Linewidth', 1.5)
%     plot(t_vec, y(i, :), '-o', 'Linewidth', 1.5)
%     yline(A(i), 'k--')
%     xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', 13)
%     ylabel('Population size', 'Interpreter', 'Latex', 'Fontsize', 13)
%     title(['Reef ', num2str(i)], 'Interpreter', 'Latex', 'Fontsize', 14) 
%     if i == 1
%         ylim([0 160])
%     end
%     if i == 2
%         legend('Coral cover (m^2)', 'No. of starfish', ...
%                'Carrying capacity of coral', 'Location', 'NorthEast')
%         ylim([0 100]) 
%     end
% end

% Predator-prey phase plane
% figure(2), clf, hold on
% for i = 1:num_reefs
%     subplot(1, num_reefs, i), hold on
%     plot(x(i, :), y(i, :), 'Linewidth', 1)
%     xlabel('Coral cover ($m^2$)', 'Interpreter', 'Latex', 'Fontsize', 13)
%     ylabel('Number of starfish', 'Interpreter', 'Latex', 'Fontsize', 13)
%     title(['Reef ', num2str(i)], 'Interpreter', 'Latex', 'Fontsize', 14)
%     if i == 1
%         xlim([0 120])
%     end
%     if i == 2
%         xlim([0 80]) 
%     end
% end

%% Plotting - Control Effort
% Predator and prey populations over time
figure(3), clf, hold on
for i = 1:num_reefs
    subplot(1, num_reefs, i), hold on
    plot(t_vec, x(i, :), '-o', 'Linewidth', 1.5)
    plot(t_vec, y(i, :), '-o', 'Linewidth', 1.5)
    yline(A(i), 'k--')
    xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', 13)
    ylabel('Population size', 'Interpreter', 'Latex', 'Fontsize', 13)
    title(['Reef ', num2str(i)], 'Interpreter', 'Latex', 'Fontsize', 14) 
    if i == 1
        ylim([0 120])
    end
    if i == 2
        legend('Coral cover (m^2)', 'No. of starfish', ...
               'Carrying capacity of coral', 'Location', 'NorthEast')
        ylim([0 100]) 
    end
end