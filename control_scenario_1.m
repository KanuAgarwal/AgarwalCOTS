% SCENARIO 1: Control at initiation box for 50 years
clear all

% MODEL ===================================================================
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

% PARAMETERS --------------------------------------------------------------
% How long do we want to run the simulation for
t_end = 100;                     % time in years

% Get the number of reefs from the lat, long data
num_reefs = length(lat);

% Parameters struct: store all parameters in a struct
% Constant parameter values for all reefs
% Estimated by Morello et al. (2014)
params.p_tilde = 0.258;         % effect of fast-growing coral on COTS     
params.M_cots = 2.56;           % natural mortality of COTS 
params.p_1_f = 0.129/2500;      % effect of COTS on fast-growing coral

% Known or arbitrarily chosen by Morello et al. (2014)
params.r_f = 0.5;               % intrinsic growth rate of fast-growing coral
params.K_f = 1;                 % carrying capacity of fast-growing coral
params.p_2_f = 10/2500;         % effect of COTS on fast-growing coral

% Known or arbitrarily chosen by me
params.r_c = 0.1;               % coral larvae reproduction rate
params.r_s = 5000;              % starfish larvae reproduction rate

% Connectivity matrices from Bode et al. (2012)
params.omega_c = psurv_d02_1122_P7;             % coral
params.omega_s = psurv_d02_1122_P7;             % starfish

% Latitude and longitude for starfish larval calculation
params.lon = lon;
params.lat = lat;

% Using metapopulation model equation for larval dispersal
dispersal_eq = 1;


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
        initial_state.N_0_2(i) = 50;
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

% Print the total number of reefs and 'budget' for control scenario
num_reefs_control_s1 = nnz(control_effort_s1(:, 1))
budget_s1 = sum(control_effort_s1, 'all')


% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec, C_y_f, N_y_2, N_y_1, N_y_0, ~] ...
    = simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s1, dispersal_eq);

% Calculate coral cover and cots over time
coral_s1 = sum(C_y_f, 1);
starfish_age2_s1 = sum(N_y_2, 1);
starfish_age1_s1 = sum(N_y_1, 1);
starfish_age0_s1 = sum(N_y_0, 1);

% Count the number of reefs with less than 1% coral 
coral_compare_all(1) = sum(C_y_f(:, end) < 0.01);

% Count the number of reefs with between 1% and 5% coral
coral_compare_all(2) = sum(C_y_f(:, end) >= 0.01 & C_y_f(:, end) < 0.05);

% Count the number of reefs with between 5% and 30% coral
coral_compare_all(3) = sum(C_y_f(:, end) >= 0.05 & C_y_f(:, end) < 0.3);

% Count the number of reefs with more than 30%
coral_compare_all(4) = sum(C_y_f(:, end) >= 0.3);

% PLOTS ===================================================================
% Define colours for plotting
green = [0.4660 0.6740 0.1880];
red = [0.8500 0.3250 0.0980];
orange = [0.9290 0.6940 0.1250];

% Fontsizes for plotting
axis_FS = 14;
title_FS = 15;
legend_FS = 13;
ticks_FS = 12;

% Define viridis and magma palette colours 
viridis_palette_4 = [253, 231, 37; 53, 183, 121; 49, 104, 142; 68, 1, 84];
viridis_palette_4 = viridis_palette_4/255;


% % Fast-growing coral heatmap ----------------------------------------------
% figure(1), clf, hold on
% imagesc(C_y_f)
% title('Fast-growing coral cover')
% xlabel('Time (years)')
% ylabel('Reef')
% colorbar
% xlim([1 t_end+1])
% ylim([0 num_reefs])

% % Age 2+ COTS heatmap -----------------------------------------------------
% figure(2), clf, hold on
% imagesc(N_y_2)
% title('Age 2+ starfish')
% xlabel('Time (years)')
% ylabel('Reef')
% colorbar
% xlim([1 t_end+1])
% ylim([0 num_reefs])

% % Total coral cover over time ---------------------------------------------
% figure(3), clf, hold on, grid on
% plot(t_vec, coral_s1, 'Linewidth', 2)
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('Coral cover (\% of reef area)', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title('Total coral cover on GBR', 'Interpreter', 'Latex', ...
%     'Fontsize', title_FS)

% % Total starfish over time ------------------------------------------------
% figure(4), clf, hold on, grid on
% plot(t_vec, starfish_age2_s1, 'Linewidth', 2)
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 2+ (adult) starfish', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title('Total adult starfish population on GBR', 'Interpreter', 'Latex', ...
%     'Fontsize', title_FS)

% % All starfish ------------------------------------------------------------
% figure(5), clf, hold on, grid on
% plot(t_vec, starfish_age2_s1, 'Linewidth', 2)
% plot(t_vec, starfish_age1_s1, 'Linewidth', 2)
% plot(t_vec, starfish_age0_s1, 'Linewidth', 2)
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of starfish', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title('Total starfish population on GBR', 'Interpreter', 'Latex', ...
%     'Fontsize', title_FS)
% legend('Age 2+ (adult)', 'Age 1 (juvenile)', 'Age 0 (larvae)', ...
%     'Interpreter', 'Latex', 'Fontsize', legend_FS, 'Location', 'NorthEastOutside')

% % Age 1 starfish ----------------------------------------------------------
% figure(5), clf, hold on, grid on
% plot(t_vec, starfish_age1_s1, 'Linewidth', 2)
% xlabel('Time (years)')
% ylabel('Total starfish')
% title('Total age 1 starfish on GBR over time')

% % Age 0 starfish ----------------------------------------------------------
% figure(6), clf, hold on, grid on
% plot(t_vec, starfish_age0_s1, 'Linewidth', 2)
% xlabel('Time (years)')
% ylabel('Total starfish')
% title('Total age 0 starfish on GBR over time')

% % Outbreak initiation on GBR ----------------------------------------------
% figure(7), clf, hold on, box on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1]);
% % Plot reef locations by color depending on initial cots numbers
% for i = 1:num_reefs
%     if initial_state.N_0_2(i) > 0
%         plot(lat(i), lon(i), 'b.', 'Markersize', 10)
%     else
%         plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', 0.7.*ones(1,3))
%     end
% end
% % Focus the figure on GBR and QLD
% xlim([140, 155])
% ylim([-26, -8])
% title('Starfish outbreak location in initiaion box')

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
figure(38), clf, hold on, box on
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
xlim([140, 155])
ylim([-26, -8])
% Plot reef locations by colour based on coral presence
for i = 1:length(lat)
    if C_y_f(i, end) >= 0.3
        p1 = plot(lat(i), lon(i), '.', 'Markersize', 12, 'Color', viridis_palette_4(1, :));
    elseif C_y_f(i, end) >= 0.05 && C_y_f(i, end) < 0.3
        p2 = plot(lat(i), lon(i), '.', 'Markersize', 12, 'Color', viridis_palette_4(2, :));
    elseif C_y_f(i, end) >= 0.01 && C_y_f(i, end) < 0.05
        p3 = plot(lat(i), lon(i), '.', 'Markersize', 12, 'Color', viridis_palette_4(3, :));
    else
        p4 = plot(lat(i), lon(i), '.', 'Markersize', 12, 'Color', viridis_palette_4(4, :));
    end
end
% Add labels
set(gca, 'FontSize', ticks_FS);
title(['Coral cover after ', num2str(t_end), ' years with 100\% effort at 168 reefs'], ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
[h, icons] = legend([p1 p2 p3 p4], '$>$30\% coral cover', '5$-$30\% coral cover', '1$-$5\% coral cover', ...
    '$<$1\% coral cover', 'Interpreter', 'Latex', 'Fontsize', legend_FS);
icons = findobj(icons, 'Type', 'line');
icons = findobj(icons, 'Marker', 'none', '-xor');
set(icons, 'MarkerSize', 25)

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

% % Starfish population on GBR ----------------------------------------------
% figure(9), clf, hold on, box on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% xlim([140, 155])
% ylim([-26, -8])
% % Plot reef locations by colour based on starfish presence
% for i = 1:length(lat)
%     if N_y_2(i, end) > 10
%         p4 = plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', red);
%     elseif N_y_2(i, end) < 0.5
%         p5 = plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', green);
%     else
%         p6 = plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', orange);
%     end
% end
% % Add labels
% set(gca, 'FontSize', ticks_FS, 'BoxStyle', 'Full');
% title(['Adult starfish population on GBR after ', num2str(t_end), ' years'], ...
%     'Interpreter', 'Latex', 'Fontsize', title_FS)
% legend([p5 p6 p4], 'No starfish', '10 or less starfish', 'More than 10 starfish', ...
%     'Interpreter', 'Latex', 'Fontsize', legend_FS)

% % Reefs being controlled on GBR -------------------------------------------
% figure(10), clf, hold on, box on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% % Plot reef locations by color depending on initial cots numbers
% for i = 1:num_reefs
%     if control_effort_s1(i, 1) > 0
%         plot(lat(i), lon(i), 'r.', 'Markersize', 10)
%     else
%         plot(lat(i), lon(i), '.', 'Markersize', 10, 'Color', 0.7.*ones(1,3))
%     end
% end
% % Focus the figure on GBR and QLD
% xlim([140, 155])
% ylim([-26, -8])
% title('Starfish control locations on the GBR')