%% SETUP
clear all

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


%% SCENARIO 1: Control at initiation box for 50 years

% MODEL ===================================================================
% PARAMETERS --------------------------------------------------------------
% How long do we want to run the simulation for
t_end = 50;                     % time in years

% Get the number of reefs from the lat, long data
num_reefs = length(lat);

% Parameters struct: store all parameters in a struct
% Constant parameter values for all reefs
% Estimated by Morello et al. (2014)
params.p_tilde = 0.258;         % effect of fast-growing coral on COTS     
params.M_cots = 2.56;           % natural mortality of COTS 
params.p_1_f = 0.129;           % effect of COTS on fast-growing coral

% Known or arbitrarily chosen by Morello et al. (2014)
params.r_f = 0.5;               % intrinsic growth rate of fast-growing coral
params.K_f = 1;                 % carrying capacity of fast-growing coral
params.p_2_f = 10;              % effect of COTS on fast-growing coral

% Known or arbitrarily chosen by me
params.r_c = 0.1;               % coral larvae reproduction rate
params.r_s = 5000;              % starfish larvae reproduction rate

% Connectivity matrices from Bode et al. (2012)
params.omega_c = psurv_d02_1122_P7;             % coral
params.omega_s = psurv_d02_1122_P7;             % starfish


% INITIAL SYSTEM STATE ----------------------------------------------------
% CORAL
% Percentage of fast-growing coral = 80% everywhere
initial_state.C_0_f = 0.8 * params.K_f * ones(num_reefs, 1);

% STARFISH
% Number of COTS aged 2+
initial_state.N_0_2 = zeros(num_reefs, 1);

% Look for reefs within the initiation box, and put some starfish there
for i = 1:num_reefs
    if (lon(i) > -17 && lon(i) < -14.75) && (lat(i) > 145 && lat(i) < 147)
        initial_state.N_0_2(i) = 100;
    end
end

% Initialise age 1 and age 0 COTS based on Morello initial conditions
initial_state.N_0_1 = initial_state.N_0_2 * exp(params.M_cots);
initial_state.N_0_0 = initial_state.N_0_2 * exp(2*params.M_cots);


% CONTROL EFFORT ----------------------------------------------------------
% Cull all the starfish only in the initiation box for all 50 years
control_effort_s1 = zeros(num_reefs, t_end);
for i = 1:num_reefs
    if (lon(i) > -17 && lon(i) < -14.75) && (lat(i) > 145 && lat(i) < 147)
        control_effort_s1(i, :) = 1;
    end
end


% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec, C_y_f, N_y_2, N_y_1, N_y_0] = simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s1);

% Calculate coral cover and cots over time
coral_s1 = sum(C_y_f, 1);
starfish_age2_s1 = sum(N_y_2, 1);
starfish_age1_s1 = sum(N_y_1, 1);
starfish_larvae_s1 = sum(N_y_0, 1);


% PLOTS ===================================================================
% Fast-growing coral heatmap ----------------------------------------------
figure(1), clf, hold on
imagesc(C_y_f)
title('Fast-growing coral cover')
xlabel('Time (years)')
ylabel('Reef')
colorbar
xlim([1 t_end+1])
ylim([0 num_reefs])

% Age 2+ COTS heatmap -----------------------------------------------------
figure(2), clf, hold on
imagesc(N_y_2)
title('Age 2+ starfish')
xlabel('Time (years)')
ylabel('Reef')
colorbar
xlim([1 t_end+1])
ylim([0 num_reefs])

% Total coral cover over time ---------------------------------------------
figure(3), clf, hold on, grid on
plot(t_vec, coral_s1, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total coral cover')
title('Total coral cover on GBR over time')

% Total starfish over time ------------------------------------------------
figure(4), clf, hold on, grid on
plot(t_vec, starfish_age2_s1, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total starfish')
title('Total age 2+ starfish on GBR over time')

figure(5), clf, hold on, grid on
plot(t_vec, starfish_age1_s1, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total starfish')
title('Total age 1 starfish on GBR over time')

figure(6), clf, hold on, grid on
plot(t_vec, starfish_larvae_s1, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total starfish')
title('Total age 0 starfish on GBR over time')


% Outbreak initiation on GBR ----------------------------------------------
figure(7), clf, hold on, box on
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
% Plot reef locations by color depending on initial cots numbers
for i = 1:num_reefs
    if initial_state.N_0_2(i) > 0
        plot(lat(i), lon(i), 'b.', 'Markersize', 10)
    else
        plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', 0.7.*ones(1,3))
    end
end
% Focus the figure on GBR and QLD
xlim([140, 155])
ylim([-26, -8])
title('Starfish outbreak location at Lizard Island')

% % GIF: Coral cover on GBR -------------------------------------------------
% h1 = figure(8); clf, hold on, box on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
% xlim([140, 155])
% ylim([-26, -8])
% % Plot reef locations by colour based on coral presence
% for t = 1:length(t_vec)
%     % Plot
%     for i = 1:length(lat)
%         if C_y_f(i, t) > 0.9
%             p(i) = plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.4660 0.6740 0.1880]);
%         elseif C_y_f(i, t) < 0.01
%             p(i) = plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.8500 0.3250 0.0980]);
%         else
%             p(i) = plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.9290 0.6940 0.1250]);
%         end
%     end
%     % Update title
%     title(['Coral cover on GBR after ', num2str(t-1), ' years'])
%     
%     % Get current plot as image
%     frame = getframe(h1);
%     im = frame2im(frame);
%     [im_ind, cm] = rgb2ind(im, 256);
% 
%     % Write to the gif file
%     if t == 1
%     elseif t == 2
%         imwrite(im_ind, cm, 'Coral.gif', 'gif', 'Loopcount', inf);
%     else
%         imwrite(im_ind, cm, 'Coral.gif', 'gif', 'WriteMode', ...
%                 'append', 'DelayTime', 0.1)
%     end
%     
%     % Draw and delete the spheres, ready for the next frame
%     drawnow
%     if t < length(t_vec)
%         delete(p)
%     end
% end

% Coral cover on GBR ------------------------------------------------------
figure(8), clf, hold on, box on
title(['Coral cover on GBR after ', num2str(t_end), ' years'])
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
xlim([140, 155])
ylim([-26, -8])
% Plot reef locations by colour based on coral presence
for i = 1:length(lat)
    if C_y_f(i, end) > 0.8
        plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.4660 0.6740 0.1880])
    elseif C_y_f(i, end) < 0.01
        plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.8500 0.3250 0.0980])
    else
        plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.9290 0.6940 0.1250])
    end
end

% % GIF: Starfish population on GBR -----------------------------------------
% h2 = figure(9); clf, hold on, box on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
% xlim([140, 155])
% ylim([-26, -8])
% % Plot reef locations by colour based on starfish presence
% for t = 1:length(t_vec)
%     % Plot
%     for i = 1:length(lat)
%         if N_y_2(i, t) > 0.5
%             p(i) = plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.8500 0.3250 0.0980]);
%         else
%             p(i) = plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.4660 0.6740 0.1880]);
%         end
%     end
%     
%     % Update title
%     title(['Starfish population on GBR after ', num2str(t-1), ' years'])
%     
%     % Get current plot as image
%     frame = getframe(h2);
%     im = frame2im(frame);
%     [im_ind, cm] = rgb2ind(im, 256);
% 
%     % Write to the gif file
%     if t == 1
%     elseif t == 2
%         imwrite(im_ind, cm, 'Starfish.gif', 'gif', 'Loopcount', inf);
%     else
%         imwrite(im_ind, cm, 'Starfish.gif', 'gif', 'WriteMode', ...
%                 'append', 'DelayTime', 0.1)
%     end
% 
%     % Draw and delete the spheres, ready for the next frame
%     drawnow
%     if t < length(t_vec)
%         delete(p)
%     end
% end

% Starfish population on GBR ----------------------------------------------
figure(9), clf, hold on, box on
title(['Starfish population on GBR after ', num2str(t_end), ' years'])
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
xlim([140, 155])
ylim([-26, -8])
% Plot reef locations by colour based on starfish presence
for i = 1:length(lat)
    if N_y_2(i, end) > 0.5
        plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.8500 0.3250 0.0980]);
    else
        plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.4660 0.6740 0.1880]);
    end
end

%% SCENARIO 2: Control at initiation box +10% buffer area for 50 years

% MODEL ===================================================================
% PARAMETERS --------------------------------------------------------------
% These are exactly the same for all scenarios, so don't need updating.

% INITIAL SYSTEM STATE ----------------------------------------------------
% These are exactly the same for all scenarios, so don't need updating.

% CONTROL EFFORT ----------------------------------------------------------
% Count number of reefs controlled in control scenario 1
num_reefs_control_s1 = nnz(control_effort_s1(:, 1))
num_reefs_control_s1 = nnz(control_effort_s1)

% Cull all starfish in the initiation box + 10% area buffer for 50 years
control_effort_s2 = zeros(num_reefs, t_end);
for i = 1:num_reefs
    if (lon(i) > -17.11 && lon(i) < -14.64) && (lat(i) > 144.75 && lat(i) < 147.25)
        control_effort_s2(i, :) = 1;
    end
end

% Count number of reefs controlled in control scenario 2
num_reefs_control_s2 = nnz(control_effort_s2(:, 1))
num_reefs_control_s2 = nnz(control_effort_s2)


% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec, C_y_f, N_y_2, N_y_1, N_y_0] = simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s2);

% Calculate coral cover and cots over time
coral_s2 = sum(C_y_f, 1);
starfish_age2_s2 = sum(N_y_2, 1);
starfish_age1_s2 = sum(N_y_1, 1);
starfish_larvae_s2 = sum(N_y_0, 1);


% PLOTS ===================================================================
% Fast-growing coral heatmap ----------------------------------------------
figure(11), clf, hold on
imagesc(C_y_f)
title('Fast-growing coral cover')
xlabel('Time (years)')
ylabel('Reef')
colorbar
xlim([1 t_end+1])
ylim([0 num_reefs])

% Age 2+ COTS heatmap -----------------------------------------------------
figure(12), clf, hold on
imagesc(N_y_2)
title('Age 2+ starfish')
xlabel('Time (years)')
ylabel('Reef')
colorbar
xlim([1 t_end+1])
ylim([0 num_reefs])

% Total coral cover over time ---------------------------------------------
figure(13), clf, hold on, grid on
plot(t_vec, coral_s2, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total coral cover')
title('Total coral cover on GBR over time')

% Total starfish over time ------------------------------------------------
figure(14), clf, hold on, grid on
plot(t_vec, starfish_age2_s2, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total starfish')
title('Total age 2+ starfish on GBR over time')

figure(15), clf, hold on, grid on
plot(t_vec, starfish_age1_s2, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total starfish')
title('Total age 1 starfish on GBR over time')

figure(16), clf, hold on, grid on
plot(t_vec, starfish_larvae_s2, 'Linewidth', 2)
xlabel('Time (years)')
ylabel('Total starfish')
title('Total age 0 starfish on GBR over time')


% Outbreak initiation on GBR ----------------------------------------------
figure(17), clf, hold on, box on
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
% Plot reef locations by color depending on initial cots numbers
for i = 1:num_reefs
    if initial_state.N_0_2(i) > 0
        plot(lat(i), lon(i), 'b.', 'Markersize', 10)
    else
        plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', 0.7.*ones(1,3))
    end
end
% Focus the figure on GBR and QLD
xlim([140, 155])
ylim([-26, -8])
title('Starfish outbreak location at Lizard Island')

% % GIF: Coral cover on GBR -------------------------------------------------
% h1 = figure(18); clf, hold on, box on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
% xlim([140, 155])
% ylim([-26, -8])
% % Plot reef locations by colour based on coral presence
% for t = 1:length(t_vec)
%     % Plot
%     for i = 1:length(lat)
%         if C_y_f(i, t) > 0.9
%             p(i) = plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.4660 0.6740 0.1880]);
%         elseif C_y_f(i, t) < 0.01
%             p(i) = plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.8500 0.3250 0.0980]);
%         else
%             p(i) = plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.9290 0.6940 0.1250]);
%         end
%     end
%     % Update title
%     title(['Coral cover on GBR after ', num2str(t-1), ' years'])
%     
%     % Get current plot as image
%     frame = getframe(h1);
%     im = frame2im(frame);
%     [im_ind, cm] = rgb2ind(im, 256);
% 
%     % Write to the gif file
%     if t == 1
%     elseif t == 2
%         imwrite(im_ind, cm, 'Coral.gif', 'gif', 'Loopcount', inf);
%     else
%         imwrite(im_ind, cm, 'Coral.gif', 'gif', 'WriteMode', ...
%                 'append', 'DelayTime', 0.1)
%     end
%     
%     % Draw and delete the spheres, ready for the next frame
%     drawnow
%     if t < length(t_vec)
%         delete(p)
%     end
% end

% Coral cover on GBR ------------------------------------------------------
figure(18), clf, hold on, box on
title(['Coral cover on GBR after ', num2str(t_end), ' years'])
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
xlim([140, 155])
ylim([-26, -8])
% Plot reef locations by colour based on coral presence
for i = 1:length(lat)
    if C_y_f(i, end) > 0.8
        plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.4660 0.6740 0.1880])
    elseif C_y_f(i, end) < 0.01
        plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.8500 0.3250 0.0980])
    else
        plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.9290 0.6940 0.1250])
    end
end

% % GIF: Starfish population on GBR -----------------------------------------
% h2 = figure(19); clf, hold on, box on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
% xlim([140, 155])
% ylim([-26, -8])
% % Plot reef locations by colour based on starfish presence
% for t = 1:length(t_vec)
%     % Plot
%     for i = 1:length(lat)
%         if N_y_2(i, t) > 0.5
%             p(i) = plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.8500 0.3250 0.0980]);
%         else
%             p(i) = plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.4660 0.6740 0.1880]);
%         end
%     end
%     
%     % Update title
%     title(['Starfish population on GBR after ', num2str(t-1), ' years'])
%     
%     % Get current plot as image
%     frame = getframe(h2);
%     im = frame2im(frame);
%     [im_ind, cm] = rgb2ind(im, 256);
% 
%     % Write to the gif file
%     if t == 1
%     elseif t == 2
%         imwrite(im_ind, cm, 'Starfish.gif', 'gif', 'Loopcount', inf);
%     else
%         imwrite(im_ind, cm, 'Starfish.gif', 'gif', 'WriteMode', ...
%                 'append', 'DelayTime', 0.1)
%     end
% 
%     % Draw and delete the spheres, ready for the next frame
%     drawnow
%     if t < length(t_vec)
%         delete(p)
%     end
% end

% Starfish population on GBR ----------------------------------------------
figure(19), clf, hold on, box on
title(['Starfish population on GBR after ', num2str(t_end), ' years'])
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
xlim([140, 155])
ylim([-26, -8])
% Plot reef locations by colour based on starfish presence
for i = 1:length(lat)
    if N_y_2(i, end) > 0.5
        plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.8500 0.3250 0.0980]);
    else
        plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', [0.4660 0.6740 0.1880]);
    end
end

