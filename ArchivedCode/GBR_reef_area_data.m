clear all

%% Reef Area CSV
% This code matches the reef areas in DataSources/ReefGazette.csv to the
% reef coordinates in IdentifyKeySources/original_centroids and keeps the
% intersection of the two datasets. It also updates the larval dispersal
% matrices in IdentifyKeySources/ConnectivityMatrices_Model_A_2002_P7
% accordingly. 


% LOAD DATA ---------------------------------------------------------------
% =========================================================================
% Load Australian coast data 
% load AustOutline            % Australian border coordinates

% =========================================================================
% Load reef coordinates 
% Load in latitude longitude data for reefs used in model
load IdentifyKeySources/original_centroids

% Rename variables to lat and lon
lon = lg;
lat = lt;
clear lg lt
reef_coord_data = [lon, lat];

% =========================================================================
% Load reef areas 
% View import options in csv and preview file
opts = detectImportOptions('DataSources/ReefGazette.csv');
preview('DataSources/ReefGazette.csv', opts)

% Only select the variables that we want to import and preview
opts.SelectedVariableNames = [15:16, 4];
preview('DataSources/ReefGazette.csv', opts)

% Load in the data with the right variables
M = readmatrix('DataSources/ReefGazette.csv', opts);

% Replace all entries that have NaN with 0
M(isnan(M)) = 0;

% Remove rows that have any zeros i.e. no area data 
reef_area_data = M(all(M, 2), :);

% =========================================================================
% Get connectivity matrices 
load IdentifyKeySources/ConnectivityMatrices_Model_A_2002_P7

% =========================================================================


% MATCH REEFS AND AREAS ---------------------------------------------------
% =========================================================================
% Plot reefs from both models
figure(10), clf, hold on
scatter(lon, lat, 15, 'red')
scatter(reef_area_data(1:end-1, 1), reef_area_data(1:end-1, 2), 8, 'blue', 'filled')
legend('Kanu model reefs', 'Reef sizes from csv file')

% =========================================================================
% Find the smallest Euclidean distance i.e. for each reef with an area,
% find the closest reef in my model
[D, I_coord] = pdist2(reef_coord_data, reef_area_data(:, 1:2), ...
    'Euclidean', 'Smallest', 1);

% Sort from smallest to largest distance
[D_sorted, I_area_sorted] = sort(D);
I_coord_sorted = I_coord(I_area_sorted);

% Get the number of unique values i.e. number of unique reefs
length(unique(I_area_sorted))
length(unique(I_coord_sorted))

% Plot both models
figure(11), clf, hold on
scatter(reef_coord_data(I_coord_sorted, 1), ...
    reef_coord_data(I_coord_sorted, 2), 15, 'red')
scatter(reef_area_data(I_area_sorted, 1), ...
    reef_area_data(I_area_sorted, 2), 8, 'blue', 'filled')
legend('Kanu model reefs', 'Reef sizes from csv file')

% =========================================================================
% Remove reefs that we couldn't find a match for from my model
unique_reefs = unique(I_coord_sorted);
reef_coord_data_v2 = reef_coord_data(unique_reefs, :);

% Remove data from connectivity matrices for reefs with no area data
omega = psurv_d02_1122_P7(unique_reefs, unique_reefs);

% =========================================================================
% Find the smallest three distances i.e. for each reef in my model, find 
% the closest 3 reefs with an area
[D_v2, I_v2_area] = pdist2(reef_area_data(:, 1:2), reef_coord_data_v2, ...
    'Euclidean', 'Smallest', 3);

% Plot both models again
figure(12), clf, hold on
scatter(reef_coord_data_v2(:, 1), ...
    reef_coord_data_v2(:, 2), 15, 'red')
scatter(reef_area_data(unique(I_v2_area(1, :)), 1), ...
    reef_area_data(unique(I_v2_area(1, :)), 2), 8, 'blue', 'filled')
legend('Kanu model reefs', 'Reef sizes from csv file')

% =========================================================================
% Create array to store reefs matches
num_matched = 0;
I_v2_area_matched = zeros(1, 1);
I_v2_coord_matched = zeros(1, 1);
num_repeats = 0;
I_v2_area_repeats = zeros(1, 1);
I_v2_coord_repeats = zeros(1, 1);

% Loop over each observation from euclidean distance algorithm
for i = 1:length(I_v2_area)
    current_reef = I_v2_area(1, i);
    % If current reef hasn't already been matched and it doesn't repeat, save it
    if ismember(current_reef, I_v2_area_matched) == 0 && ismember(current_reef, I_v2_area_repeats) == 0
        num_matched = num_matched + 1; 
        I_v2_area_matched(num_matched) = current_reef;
        I_v2_coord_matched(num_matched) = i;
    % Otherwise, if current reef has already been matched or repeated
    else
        % If the reef hasn't repeated before 
        if ismember(current_reef, I_v2_area_repeats) == 0
            % Find previous entry in array
            index_remove = find(I_v2_area_matched == current_reef);
            % Add previous entry to repeat array
            num_repeats = num_repeats + 1;
            I_v2_area_repeats(num_repeats) = I_v2_area_matched(index_remove);
            I_v2_coord_repeats(num_repeats) = I_v2_coord_matched(index_remove);
            % Remove entry from matched reefs
            I_v2_area_matched(index_remove) = [];
            I_v2_coord_matched(index_remove) = [];
            num_matched = num_matched - 1;
            % Add new entry to repeat array
            num_repeats = num_repeats + 1;
            I_v2_area_repeats(num_repeats) = current_reef;
            I_v2_coord_repeats(num_repeats) = i;
        % Otherwise if the reef has repeated before
        else 
            % Add new entry to repeat array
            num_repeats = num_repeats + 1;
            I_v2_area_repeats(num_repeats) = current_reef;
            I_v2_coord_repeats(num_repeats) = i;
        end
    end
end

% =========================================================================
% Get all the distances and indexes that repeat
D_v2_repeats = D_v2(:, I_v2_coord_repeats);
I_v2_area_repeats_all = I_v2_area(:, I_v2_coord_repeats);
I_v2_area_repeats_unique = unique(I_v2_area_repeats);

% Create new array for area data
I_v2_area_actual = I_v2_area(1, :);

% Find the closest reef for any repeats
for i = 1:length(I_v2_area_repeats_unique)
    % Find all repeats of each unique element
    current_repeats = find(I_v2_area_repeats == I_v2_area_repeats_unique(i));
    % Make a new array with all for index and distance
    D_current = D_v2_repeats(:, current_repeats);
    I_area_current = I_v2_area_repeats_all(:, current_repeats);
    I_coord_current = I_v2_coord_repeats(current_repeats);
    % Find the smallest distance in the first row i.e. the closest reef
    [~, ind_min] = min(D_current(1, :));
    [~, ind_max] = max(D_current(1, :));
    % Update the index so that the closest reef is used 
    I_v2_area_actual(I_coord_current(ind_min)) = I_area_current(1, ind_min);
    % Figure out the closest reef for the other repeats
    % If there is only one other reef
    if size(D_current, 2) == 2
        % Make sure the second closest reef isn't already used
        if ~any(I_v2_area_actual == I_area_current(2, ind_max))
            I_v2_area_actual(I_coord_current(ind_max)) = I_area_current(2, ind_max);
        % Make sure the third closest reef isn't already used
        elseif ~any(I_v2_area_actual == I_area_current(3, ind_max))
            I_v2_area_actual(I_coord_current(ind_max)) = I_area_current(3, ind_max);
        % Otherwise place a zero and delete the entry later
        else 
            I_v2_area_actual(I_coord_current(ind_max)) = 0;
        end
    % If there are two other reefs with the same 
    elseif size(D_current, 2) == 3
        % Get the subset of indexes and values
        D_subset = D_current(2:3, 1:end ~= ind_min);
        I_area_subset = I_area_current(2:3, 1:end ~= ind_min);
        I_coord_subset = I_coord_current(1:end ~= ind_min);
        % Check that the second closest reefs aren't the same
        if I_area_subset(1, 1) ~= I_area_subset(1, 2)
            % Make sure the first reef isn't already used
            if ~any(I_v2_area_actual == I_area_subset(1, 1))
                I_v2_area_actual(I_coord_current(1)) = I_area_subset(1, 1);
            else
                fprintf('help reef 1 already used')
            end
            % Make sure the second reef isn't already used
            if ~any(I_v2_area_actual == I_area_subset(1, 2))
                I_v2_area_actual(I_coord_current(2)) = I_area_subset(1, 2);
            else
                fprintf('help reef 2 already used')
            end
        else
            fprintf('help same reefs again')
        end
    end
end

% Delete any entries from indexes that are zero - these are repeat reefs
index_delete = find(I_v2_area_actual == 0); 
I_v2_area_actual(index_delete) = [];

% Delete the corresponding data from everything else i.e. connectivity 
% matrix, lat/long data
reef_coord_data_v2(index_delete, :) = []; 
omega(index_delete, :) = [];
omega(:, index_delete) = [];

% =========================================================================
% Check all matching reefs are unique 
num_reefs_actual = length(reef_coord_data_v2)
unique_reefs = length(unique(I_v2_area_actual))

% Keep only the reefs that match
reef_area_data_actual = reef_area_data(I_v2_area_actual, :);

% Plot both models again
figure(13), clf, hold on
scatter(reef_coord_data_v2(:, 1), ...
    reef_coord_data_v2(:, 2), 15, 'red')
scatter(reef_area_data(I_v2_area_actual, 1), ...
    reef_area_data(I_v2_area_actual, 2), 8, 'blue', 'filled')
legend('Kanu model reefs', 'Reef sizes from csv file')

% =========================================================================
% Variables to return
num_reefs = num_reefs_actual;
lon = reef_coord_data_v2(:, 1);
lat = reef_coord_data_v2(:, 2);
reef_area = reef_area_data_actual(:, 3);
omega = omega;


%% Reef Polygons

% % Load data ---------------------------------------------------------------
% load ReefOutline            % GBR individual reef coordinates
% load AustOutline            % Australian border coordinates
% 
% % Load in latitude longitude data for reefs used in model
% load IdentifyKeySources/original_centroids
% 
% % Rename variables to lat and lon
% lon = lg;
% lat = lt;
% clear lg lt
% 
% 
% % Create reef polygons ----------------------------------------------------
% % Create a single polygon with all reef coordinates 
% GBR_polygon = polyshape(ReefRaw(:, 1), ReefRaw(:, 2));
% 
% % Get correct reef coordinates from polygon
% reef_coordinates = GBR_polygon.Vertices;
% 
% % Add NaN to the last row so each polygon can be extracted properly
% reef_coordinates(end+1, :) = NaN;
% 
% % Initialise variables to keep track of index
% index = 1;
% polygon_index = 1;
% 
% % Make a separate polygon for each individual reef on GBR
% for i = 1:size(reef_coordinates, 1)
%     % If the current set of coordinates is NaN, we are at a new reef
%     if isnan(reef_coordinates(i, 1))
%         % Save the reef coordinates
%         reef_polygons(index) = polyshape(reef_coordinates(polygon_index:i-1, 1), ...
%             reef_coordinates(polygon_index:i-1, 2));
%         % Increment the indexes 
%         index = index + 1;
%         polygon_index = i + 1;
%     end
% end
% 
% 
% % Calculate area of each reef ---------------------------------------------
% % Initalise array to store areas in
% reef_areas = zeros(1, length(reef_polygons));
% 
% % Get the earth radius in km to help calculate area
% rKms = earthRadius('km');
% 
% % Loop over each reef polygon
% for i = 1:length(reef_polygons)
%     reef_areas(i) = areaint(reef_polygons(i).Vertices(:, 2), ...
%         reef_polygons(i).Vertices(:, 1), rKms);
% end
% 
%
% % Match reefs used in model to reef polygons ------------------------------
% % Initialise array to store matched reef areas
% model_reef_areas = zeros(1, length(lat));
% 
%
% % Loop over every reef polygon
% for i = 1:length(reef_polygons)
%     % Decide whether current polygon has a point lying inside it
%     [in, on] = inpolygon(lon, lat, reef_polygons(i).Vertices(:, 1), reef_polygons(i).Vertices(:, 2));
%     
%     num_ins = nnz(in)
%     num_ons = nnz(on)
%     
%     figure(60), clf, hold on
%     plot(reef_polygons(i))
%     xl = xlim;
%     yl = ylim;
%     scatter(lon, lat, 10, 'black', 'filled')
%     scatter(lon(in), lat(in), 10, 'green', 'filled')
%     scatter(lon(on), lat(on), 10, 'blue', 'filled')
%     xlim([xl(1) xl(2)])
%     ylim([yl(1) yl(2)])
%     title(['Reef ', num2str(i)])
%     
%     % If a point is found
%     if numel(lon(in)) == 0
%         point_index = 0;
%         num_points = 0;
%         % Find the closest reef point coords if it exists
%         for j = 1:length(lon)
%             if lon(j) >= xl(1) && lon(j) <= xl(2) && lat(j) >= yl(1) && lat(j) <= yl(2)
%                 num_points = num_points + 1;
%                 point_index(num_points) = j;
%             end
%         end
%         % Assign reef area to correct index
%         if nnz(point_index) == 1
%             model_reef_areas(point_index) = reef_areas(i);
%         elseif nnz(point_index) > 1
%             point_index
%         end
%     elseif numel(lon(in)) == 1
%         % Get the index of the reef 
%         point_index = find(in == 1);
%         % Assign the reef area to the correct index
%         model_reef_areas(point_index) = reef_areas(i);
%     elseif numel(lon(in)) > 1
%         point_index = find(in == 1)
%     end
%     
%     num_areas = nnz(model_reef_areas)
%     point_index = 0;
% end


% % PLOTS ===================================================================
% % Plot polygon of entire GBR ----------------------------------------------
% figure(61), clf, hold on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% % Focus the figure on GBR and QLD
% xlim([140, 155])
% ylim([-26, -8])
% % Plot reefs 
% plot(GBR_polygon)
% scatter(lon, lat, 10, 'black', 'filled')
%
% % Plot ploygon of GBR with individual reefs -------------------------------
% figure(2), clf, hold on
% % Plot outline of Australia
% pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% % Focus the figure on GBR and QLD
% xlim([140, 155])
% ylim([-26, -8])
% % Plot reefs
% plot(reef_polygons)
% scatter(lon, lat, 10, 'black', 'filled')


%% AIMS Benthic Layer Data

% % Load in raster layer - tif file
% [A, R] = readgeoraster('DataSources/GBR10_GBRMP_Benthic/GBR10 GBRMP Benthic.tif');
% 
% % Create datastore to split up raster layer since its too large
% ds = datastore('DataSources/GBR10_GBRMP_Benthic/GBR10 GBRMP Benthic.tif', 'Type', 'Image');
% 
% % Partition datastore 
% new_ds = partition(ds, 10, 1);
% 
% % Convert to tall array
% T = tall(ds);
% 
% % Load in as tiff and read it
% t = Tiff('DataSources/GBR10_GBRMP_Benthic/GBR10 GBRMP Benthic.tif', 'r');
% imageData = read(t);
% 
% % Use imread to display
% imread('DataSources/GBR10_GBRMP_Benthic/GBR10 GBRMP Benthic.tif')
