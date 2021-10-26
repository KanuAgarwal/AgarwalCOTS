% SCENARIO 0: No control, but simulation run with same conditions
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
% No control
control_effort_s0 = 0;


% SOLVE -------------------------------------------------------------------
% Solve using function which runs simulations
[t_vec, C_y_f, N_y_2, N_y_1, N_y_0, tau_ratio] = ...
    simulate_reefs_v2(num_reefs, t_end, params, initial_state, control_effort_s0, dispersal_eq);

% Calculate coral cover and cots over time
coral_s0 = sum(C_y_f, 1);
starfish_age2_s0 = sum(N_y_2, 1);
starfish_age1_s0 = sum(N_y_1, 1);
starfish_age0_s0 = sum(N_y_0, 1);

% Calculate coral cover and starfish population in box
[coral_box, starfish_age2_box, starfish_age1_box, starfish_age0_box] ...
    = calculate_population_box(t_end, C_y_f, N_y_2, N_y_1, N_y_0, num_reefs, lat, lon);

% Count the number of reefs with less than 1% coral 
coral_compare_all(1) = sum(C_y_f(:, end) < 0.01);

% Count the number of reefs with between 1% and 5% coral
coral_compare_all(2) = sum(C_y_f(:, end) >= 0.01 & C_y_f(:, end) < 0.05);

% Count the number of reefs with between 5% and 30% coral
coral_compare_all(3) = sum(C_y_f(:, end) >= 0.05 & C_y_f(:, end) < 0.3);

% Count the number of reefs with more than 30%
coral_compare_all(4) = sum(C_y_f(:, end) >= 0.3);

%% PLOTS ===================================================================
% Define colours for plotting
green = [0.4660 0.6740 0.1880];
red = [0.8500 0.3250 0.0980];
orange = [0.9290 0.6940 0.1250];

% Define viridis and magma palette colours 
viridis_palette_3 = [253, 231, 37; 33, 145, 140; 68, 1, 84];
viridis_palette_3 = viridis_palette_3/255;
viridis_palette_4 = [253, 231, 37; 53, 183, 121; 49, 104, 142; 68, 1, 84];
viridis_palette_4 = viridis_palette_4/255;
viridis_palette_5 = [253, 231, 3; 94, 201, 98; 33, 145, 140; 59, 82, 139; 68, 1, 84];
viridis_palette_5 = viridis_palette_5/255;
viridis_palette_6 = [253, 231, 37; 122, 209, 81; 34, 168, 132; 42, 120, 142; 65, 68, 135; 68, 1, 84];
viridis_palette_6 = viridis_palette_6/255;
magma_palette = [252, 253, 191; 252, 137, 97; 183, 55, 121; 81, 18, 124; 0, 0, 4];
magma_palette = magma_palette/255;

% Color values from cbrewer
colour_scheme = cbrewer('div', 'PRGn', 9);

% Fontsizes for plotting
axis_FS = 14;
title_FS = 15;
legend_FS = 13;
ticks_FS = 12;

% Load colormaps
% jet=colormap('jet');
parula=fake_parula();
magma=magma();
inferno=inferno();
plasma=plasma();
viridis=viridis();


%% POPULAITONS ============================================================
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

% Total coral cover over time ---------------------------------------------
figure(13), clf, hold on
plot(t_vec, coral_s0, 'Linewidth', 2, 'Color', colour_scheme(1, :))
yline(num_reefs, '--', 'Linewidth', 2, 'Color', [0.5 0.5 0.5])
set(gca, 'FontSize', ticks_FS)
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total coral cover', 'Interpreter', 'Latex', ...
    'Fontsize', axis_FS)
title('Total coral cover on GBR', 'Interpreter', 'Latex', ...
    'Fontsize', title_FS)
legend('Total coral cover', 'Carrying capacity', 'Interpreter', 'Latex', ... 
    'Fontsize', legend_FS-1, 'Location', 'SouthWest')

% Total coral cover at initiation box over time ---------------------------
figure(14), clf, hold on
plot(t_vec, coral_box, 'Linewidth', 2, 'Color', colour_scheme(1, :))
yline(nnz(initial_state.N_0_2), '--', 'Linewidth', 2, 'Color', [0.5 0.5 0.5])
set(gca, 'FontSize', ticks_FS)
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('Total coral cover', 'Interpreter', 'Latex', ...
    'Fontsize', axis_FS)
title('Total coral cover at initiation box', 'Interpreter', 'Latex', ...
    'Fontsize', title_FS)
legend('Total coral cover', 'Carrying capacity', 'Interpreter', 'Latex', ... 
    'Fontsize', legend_FS-1, 'Location', 'SouthWest')

% % Total age 2+ starfish over time -----------------------------------------
% figure(14), clf, hold on, grid on
% plot(t_vec, starfish_age2_s0, 'Linewidth', 2)
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of age 2+ (adult) starfish', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title('Total adult starfish population on GBR', 'Interpreter', 'Latex', ...
%     'Fontsize', title_FS)

% % Total age 1 starfish over time ------------------------------------------
% figure(5), clf, hold on, grid on
% plot(t_vec, starfish_age1_s0, 'Linewidth', 2)
% xlabel('Time (years)')
% ylabel('Total starfish')
% title('Total age 1 starfish on GBR over time')

% % Total age 0 starfish over time ------------------------------------------
% figure(6), clf, hold on, grid on
% plot(t_vec, starfish_age0_s0, 'Linewidth', 2)
% xlabel('Time (years)')
% ylabel('Total starfish')
% title('Total age 0 starfish on GBR over time')

% % All starfish ------------------------------------------------------------
% figure(15), clf, hold on, grid on
% plot(t_vec, starfish_age2_s0, 'Linewidth', 2)
% plot(t_vec, starfish_age1_s0, 'Linewidth', 2)
% plot(t_vec, starfish_age0_s0, 'Linewidth', 2)
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('No. of starfish', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title('Total starfish population on GBR', 'Interpreter', 'Latex', ...
%     'Fontsize', title_FS)
% legend('Age 2+ (adult)', 'Age 1 (juvenile)', 'Age 0 (larvae)', ...
%     'Interpreter', 'Latex', 'Fontsize', legend_FS, 'Location', 'NorthWest')

% Log of all starfish -----------------------------------------------------
figure(16), clf, hold on
plot(t_vec, starfish_age0_s0, 'Linewidth', 2, 'Color', colour_scheme(7, :))
plot(t_vec, starfish_age1_s0, 'Linewidth', 2, 'Color', colour_scheme(8, :))
plot(t_vec, starfish_age2_s0, 'Linewidth', 2, 'Color', colour_scheme(9, :))
set(gca, 'YScale', 'log')
set(gca, 'FontSize', ticks_FS);
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('No. of starfish', 'Interpreter', 'Latex', ...
    'Fontsize', axis_FS)
title('Total starfish population on GBR', 'Interpreter', 'Latex', ...
    'Fontsize', title_FS)
legend('Larvae (age 0)', 'Juvenile (age 1)', 'Adult (age 2+)', ...
    'Interpreter', 'Latex', 'Fontsize', legend_FS-1, 'Location', 'NorthWest')

% Log of all starfish in box ----------------------------------------------
figure(15), clf, hold on
plot(t_vec, starfish_age0_box, 'Linewidth', 2, 'Color', colour_scheme(7, :))
plot(t_vec, starfish_age1_box, 'Linewidth', 2, 'Color', colour_scheme(8, :))
plot(t_vec, starfish_age2_box, 'Linewidth', 2, 'Color', colour_scheme(9, :))
set(gca, 'YScale', 'log')
set(gca, 'FontSize', ticks_FS);
xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
ylabel('No. of starfish', 'Interpreter', 'Latex', ...
    'Fontsize', axis_FS)
title('Total starfish population at initiation box', 'Interpreter', 'Latex', ...
    'Fontsize', title_FS)
legend('Larvae (age 0)', 'Juvenile (age 1)', 'Adult (age 2+)', ...
    'Interpreter', 'Latex', 'Fontsize', legend_FS-1, 'Location', 'NorthWest')


% % Subplots ----------------------------------------------------------------
% figure(17), clf, hold on
% 
% % Coral over time
% subplot(2, 2, 1), hold on 
% plot(t_vec, coral_s0, 'Linewidth', 2, 'Color', [0.4940 0.1840 0.5560])
% set(gca, 'FontSize', ticks_FS);
% ylabel('Total coral cover', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title('Total coral cover on GBR', 'Interpreter', 'Latex', ...
%     'Fontsize', title_FS)
% 
% % Coral in box over time
% subplot(2, 2, 2), hold on
% plot(t_vec, coral_box, 'Linewidth', 2, 'Color', [0.4940 0.1840 0.5560])
% set(gca, 'FontSize', ticks_FS);
% title('Total coral cover at initiation box', 'Interpreter', 'Latex', ...
%     'Fontsize', title_FS)
% 
% % All starfish over time
% subplot(2, 2, 3), hold on
% plot(t_vec, starfish_age2_s0, 'Linewidth', 2)
% plot(t_vec, starfish_age1_s0, 'Linewidth', 2)
% plot(t_vec, starfish_age0_s0, 'Linewidth', 2)
% set(gca, 'YScale', 'log')
% set(gca, 'YGrid', 'on')
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% ylabel('Total no. of starfish', 'Interpreter', 'Latex', ...
%     'Fontsize', axis_FS)
% title('Total starfish population on GBR', 'Interpreter', 'Latex', ...
%     'Fontsize', title_FS)
% 
% % All starfish in box over time
% subplot(2, 2, 4), hold on 
% plot(t_vec, starfish_age2_box, 'Linewidth', 2)
% plot(t_vec, starfish_age1_box, 'Linewidth', 2)
% plot(t_vec, starfish_age0_box, 'Linewidth', 2)
% set(gca, 'YScale', 'log')
% set(gca, 'YGrid', 'on')
% set(gca, 'FontSize', ticks_FS);
% xlabel('Time (years)', 'Interpreter', 'Latex', 'Fontsize', axis_FS)
% title('Total starfish population at initiation box', 'Interpreter', 'Latex', ...
%     'Fontsize', title_FS)
% legend('Age 2+ (adult)', 'Age 1 (juvenile)', 'Age 0 (larvae)', ...
%     'Interpreter', 'Latex', 'Fontsize', legend_FS-1, 'Location', 'NorthWest')


%% MAPS ===================================================================
% Outbreak initiation on GBR ----------------------------------------------
figure(17), clf, hold on, box on
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% Plot reef locations by color depending on initial cots numbers
for i = 1:num_reefs
    if initial_state.N_0_2(i) > 0
        pr = plot(lat(i), lon(i), '.', 'Markersize', 12, 'Color', magma_palette(4, :));
    else
        pg = plot(lat(i), lon(i), '.', 'Markersize', 12, 'Color', magma_palette(5, :));
    end
end
% Focus the figure on GBR and QLD
xlim([140, 155])
ylim([-26, -8])
% Add labels
set(gca, 'FontSize', ticks_FS);
title('Initial adult starfish population at initiation box', 'Interpreter', 'Latex', ...
    'Fontsize', title_FS)
[h, icons] = legend([pr pg], '$<$ 10$^2$ starfish', '0 starfish', 'Interpreter', 'Latex', 'Fontsize', legend_FS);
icons = findobj(icons, 'Type', 'line');
icons = findobj(icons, 'Marker', 'none', '-xor');
set(icons, 'MarkerSize', 25)

% % Outbreak initiation on GBR (colormap) -----------------------------------
% figure(17), clf, hold on, box on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% % Plot reef locations by color depending on initial cots numbers
% scatter(lat, lon, 12, N_y_2(:, 1), 'filled')
% colorbar
% colormap(magma)
% caxis([0 1.8812e+08])
% % Focus the figure on GBR and QLD
% xlim([140, 155])
% ylim([-26, -8])
% % Add labels
% set(gca, 'FontSize', ticks_FS, 'BoxStyle', 'Full');
% title('Starfish outbreak locations in initiation box', 'Interpreter', 'Latex', ...
%     'Fontsize', title_FS)

% % GIF: Coral cover on GBR -------------------------------------------------
% h1 = figure(19); clf, hold on, box on
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
title(['Coral cover after ', num2str(t_end), ' years with no control'], ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
[h, icons] = legend([p1 p2 p3 p4], '$>$30\% coral cover', '5$-$30\% coral cover', '1$-$5\% coral cover', ...
    '$<$1\% coral cover', 'Interpreter', 'Latex', 'Fontsize', legend_FS);
icons = findobj(icons, 'Type', 'line');
icons = findobj(icons, 'Marker', 'none', '-xor');
set(icons, 'MarkerSize', 25)


% Coral cover on GBR (colormap) -------------------------------------------
figure(19), clf, hold on, box on
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% Plot reef locations by color depending on initial cots numbers
scatter(lat, lon, 12, C_y_f(:, end), 'filled')
colorbar
colormap(viridis)
% Focus the figure on GBR and QLD
xlim([140, 155])
ylim([-26, -8])
% Add labels
set(gca, 'FontSize', ticks_FS, 'BoxStyle', 'Full');
title(['Coral cover after ', num2str(t_end), ' years with no control'], ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)

% % GIF: Starfish population on GBR -----------------------------------------
% h2 = figure(20); clf, hold on, box on
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
figure(20), clf, hold on, box on
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
xlim([140, 155])
ylim([-26, -8])
% Plot reef locations by colour based on starfish presence
for i = 1:length(lat)
    if N_y_2(i, end) >= 10^6
        p7 = plot(lat(i), lon(i), '.', 'Markersize', 12, 'Color', magma_palette(1, :));
    elseif N_y_2(i, end) >= 10^4 && N_y_2(i, end) < 10^6
        p8 = plot(lat(i), lon(i), '.', 'Markersize', 12, 'Color', magma_palette(2, :));
    elseif N_y_2(i, end) >= 10^2 && N_y_2(i, end) < 10^4
        p9 = plot(lat(i), lon(i), '.', 'Markersize', 12, 'Color', magma_palette(3, :));
    elseif N_y_2(i, end) >= 0.5 && N_y_2(i, end) < 10^2
        p10 = plot(lat(i), lon(i), '.', 'Markersize', 12, 'Color', magma_palette(4, :));
    else
        p11 = plot(lat(i), lon(i), '.', 'Markersize', 12, 'Color', magma_palette(5, :));
    end
end
% Add labels
set(gca, 'FontSize', ticks_FS, 'BoxStyle', 'Full');
title(['Adult starfish population after ', num2str(t_end), ' years with no control'], ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)
[h, icons] = legend([p7 p8 p9 p10 p11], '$>$ 10$^6$ starfish', '10$^4-$10$^6$ starfish', '10$^2-$10$^4$ starfish', ...
    '$<$ 10$^2$ starfish', '0 starfish', 'Interpreter', 'Latex', 'Fontsize', legend_FS);
icons = findobj(icons, 'Type', 'line');
icons = findobj(icons, 'Marker', 'none', '-xor');
set(icons, 'MarkerSize', 25)

% Starfish population on GBR (colormap) -----------------------------------
figure(21), clf, hold on, box on
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% Plot reef locations by color depending on initial cots numbers
scatter(lat, lon, 12, N_y_2(:, end), 'filled')
colorbar
colormap(magma)
% Focus the figure on GBR and QLD
xlim([140, 155])
ylim([-26, -8])
% Add labels
set(gca, 'FontSize', ticks_FS, 'BoxStyle', 'Full');
title(['Adult starfish population after ', num2str(t_end), ' years with no control'], ...
    'Interpreter', 'Latex', 'Fontsize', title_FS)