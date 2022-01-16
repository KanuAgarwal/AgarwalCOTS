clear all

%% Reef Polygons

% Load data ---------------------------------------------------------------
load ReefOutline            % GBR individual reef coordinates
load AustOutline            % Australian border coordinates

% Load in latitude longitude data for reefs used in model
load IdentifyKeySources/original_centroids

% Rename variables to lat and lon
lon = lg;
lat = lt;
clear lg lt


% Create reef polygons ----------------------------------------------------
% Create a single polygon with all reef coordinates 
GBR_polygon = polyshape(ReefRaw(:, 1), ReefRaw(:, 2));

% Get correct reef coordinates from polygon
reef_coordinates = GBR_polygon.Vertices;

% Add NaN to the last row so each polygon can be extracted properly
reef_coordinates(end+1, :) = NaN;

% Initialise variables to keep track of index
index = 1;
polygon_index = 1;

% Make a separate polygon for each individual reef on GBR
for i = 1:size(reef_coordinates, 1)
    % If the current set of coordinates is NaN, we are at a new reef
    if isnan(reef_coordinates(i, 1))
        % Save the reef coordinates
        reef_polygons(index) = polyshape(reef_coordinates(polygon_index:i-1, 1), ...
            reef_coordinates(polygon_index:i-1, 2));
        % Increment the indexes 
        index = index + 1;
        polygon_index = i + 1;
    end
end


% Calculate area of each reef ---------------------------------------------
% Initalise array to store areas in
reef_areas = zeros(1, length(reef_polygons));

% Get the earth radius in km to help calculate area
rKms = earthRadius('km');

% Loop over each reef polygon
for i = 1:length(reef_polygons)
    reef_areas(i) = areaint(reef_polygons(i).Vertices(:, 2), ...
        reef_polygons(i).Vertices(:, 1), rKms);
end


% Match reefs used in model to reef polygons ------------------------------
% Initialise array to store matched reef areas
model_reef_areas = zeros(1, length(lat));

% Loop over every reef polygon
for i = 1:length(reef_polygons)
    % Decide whether current polygon has a point lying inside it
    in = inpolygon(lon, lat, reef_polygons(i).Vertices(:, 1), reef_polygons(i).Vertices(:, 2));
    
    % If a point is found
    if numel(lon(in)) == 1
        % Get the index of the reef 
        point_index = find(in == 1);
        % Assign the reef area to the correct index
        model_reef_areas(point_index) = reef_areas(i);
    elseif numel(lon(in)) > 1
        point_index = find(in == 1)
    end
end


% PLOTS ===================================================================
% Plot polygon of entire GBR ----------------------------------------------
figure(1), clf, hold on
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% Focus the figure on GBR and QLD
xlim([140, 155])
ylim([-26, -8])
% Plot reefs 
plot(GBR_polygon)

% Plot ploygon of GBR with individual reefs -------------------------------
figure(2), clf, hold on
% Plot outline of Australia
pt = patch(Outline(:, 1), Outline(:, 2), [1 1 1], 'FaceColor', [0.8 0.8 0.8]);
% Focus the figure on GBR and QLD
xlim([140, 155])
ylim([-26, -8])
% Plot reefs
plot(reef_polygons)
scatter(lat, lon, 50)


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
