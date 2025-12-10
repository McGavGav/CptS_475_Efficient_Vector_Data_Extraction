// Load Washington counties
var counties = ee.FeatureCollection('TIGER/2018/Counties');
var waCounties = counties.filter(ee.Filter.eq('STATEFP', '53'));
Map.centerObject(waCounties, 6);
Map.addLayer(waCounties, {}, 'WA Counties');

// STRATEGY 1: Spatial Chunking/Tiling
// Divides large polygons into manageable grid cells

function createSpatialGrid(geometry, cellSize) {
  // cellSize in degrees (e.g., 0.1 = ~11km at this latitude)
  var bounds = geometry.bounds();
  var coords = ee.List(bounds.coordinates().get(0));
  var xMin = ee.List(coords.get(0)).get(0);
  var yMin = ee.List(coords.get(0)).get(1);
  var xMax = ee.List(coords.get(2)).get(0);
  var yMax = ee.List(coords.get(2)).get(1);
  
  var xSteps = ee.Number(xMax).subtract(xMin).divide(cellSize).ceil();
  var ySteps = ee.Number(yMax).subtract(yMin).divide(cellSize).ceil();
  
  var cells = ee.FeatureCollection(
    ee.List.sequence(0, xSteps.subtract(1)).map(function(x) {
      return ee.List.sequence(0, ySteps.subtract(1)).map(function(y) {
        var xStart = ee.Number(xMin).add(ee.Number(x).multiply(cellSize));
        var yStart = ee.Number(yMin).add(ee.Number(y).multiply(cellSize));
        var cell = ee.Geometry.Rectangle([
          xStart, yStart,
          xStart.add(cellSize), yStart.add(cellSize)
        ]);
        return ee.Feature(cell.intersection(geometry, 1), {
          'grid_x': x,
          'grid_y': y
        });
      });
    }).flatten()
  );
  
  return cells.filter(ee.Filter.neq('area', 0));
}

// STRATEGY 2: Hierarchical Processing
// Process by administrative levels (state -> county -> smaller units)

function processCountyInChunks(county, envFunction) {
  var grid = createSpatialGrid(county.geometry(), 0.05); // ~5.5km cells
  
  var results = grid.map(function(cell) {
    return envFunction(cell);
  });
  
  return results;
}

// STRATEGY 3: Adaptive Scale Reduction
// Automatically adjusts scale based on polygon size

function getAdaptiveScale(geometry) {
  var area = geometry.area(1); // in square meters
  // Larger areas get coarser resolution to avoid timeout
  var scale = ee.Number(area).sqrt().divide(1000).max(30).min(500);
  return scale;
}

// PART 2: ENVIRONMENTAL DATA EXTRACTION

// Function 1: NDVI (Green Space Indicator) - Your existing work
function computeNDVI(feature, startDate, endDate) {
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(feature.geometry())
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10))
    .select(['B4', 'B8'])
    .median();
  
  var ndvi = s2.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndviClipped = ndvi.clip(feature.geometry());
  
  return ndviClipped;
}

// Function 2: Air Quality (NO2 - Nitrogen Dioxide)
function computeNO2(feature, startDate, endDate) {
  var no2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2')
    .filterBounds(feature.geometry())
    .filterDate(startDate, endDate)
    .select('tropospheric_NO2_column_number_density')
    .mean()
    .clip(feature.geometry());
  
  return no2;
}

// Function 3: Temperature (Land Surface Temperature)
function computeLST(feature, startDate, endDate) {
  var lst = ee.ImageCollection('MODIS/061/MOD11A2')
    .filterBounds(feature.geometry())
    .filterDate(startDate, endDate)
    .select('LST_Day_1km')
    .mean()
    .multiply(0.02) // Scale factor
    .subtract(273.15) // Convert Kelvin to Celsius
    .clip(feature.geometry());
  
  return lst;
}

// Function 4: PM2.5 (Particulate Matter)
// OPTION 1: NASA GEOS-CF
function computePM25_GEOSCF(feature, startDate, endDate) {
  var pm25 = ee.ImageCollection('NASA/GEOS-CF/v1/rpl/tavg1hr')
    .filterBounds(feature.geometry())
    .filterDate(startDate, endDate)
    .select('PM25_RH35_GCC') // PM2.5 at 35% relative humidity
    .mean()
    .clip(feature.geometry());
  
  return pm25;
}
// OPTION 2: Washington University Global PM2.5 (Annual/Monthly averages)
function computePM25_WashU(feature, year) {
  var pm25 = ee.ImageCollection('projects/sat-io/open-datasets/GLOBAL_SATELLITE_PM25')
    .filter(ee.Filter.eq('year', year))
    .first()
    .clip(feature.geometry());
  
  return pm25;
}

// OPTION 3: GHAP Daily PM2.5 (2017-2022, 1km resolution)
function computePM25_GHAP(feature, startDate, endDate) {
  var pm25 = ee.ImageCollection('projects/sat-io/open-datasets/GHAP/GHAP_D1K_PM25')
    .filterBounds(feature.geometry())
    .filterDate(startDate, endDate)
    .mean()
    .clip(feature.geometry());
  
  return pm25;
}

// Default function - uses GEOS-CF for flexibility with date ranges
function computePM25(feature, startDate, endDate) {
  return computePM25_GEOSCF(feature, startDate, endDate);
}

// Function to create points along line for more detailed sampling
function createPointsAlongLine(lineGeometry, intervalMeters) {
  var length = lineGeometry.length();
  var distances = ee.List.sequence(0, length, intervalMeters);
  
  var points = distances.map(function(dist) {
    return lineGeometry.cutLines([dist]).geometries().get(0);
  });
  
  return ee.FeatureCollection(points.map(function(geom) {
    return ee.Feature(geom);
  }));
}

// Advanced exposure calculation with spatial resolution
function calculateDetailedExposure(trip, envLayer, samplingInterval) {
  var points = createPointsAlongLine(trip.geometry(), samplingInterval);
  
  var sampledPoints = points.map(function(point) {
    var value = envLayer.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: point.geometry().buffer(25),
      scale: 30,
      maxPixels: 1e9
    });
    
    return point.set('exposure', value.values().get(0));
  });
  
  var avgExposure = sampledPoints.aggregate_mean('exposure');
  var maxExposure = sampledPoints.aggregate_max('exposure');
  var minExposure = sampledPoints.aggregate_min('exposure');
  
  return trip.set({
    'mean_exposure': avgExposure,
    'max_exposure': maxExposure,
    'min_exposure': minExposure,
    'sampling_points': sampledPoints.size()
  });
}

// PART 4: INTEGRATED WORKFLOW - ENVIRONMENTAL LAYERS 
// Set date range for environmental data
var startDate = '2023-05-01';
var endDate = '2023-06-30';

// Create environmental layers for entire Washington State
var waGeometry = waCounties.geometry();

print('Computing environmental layers...');
print('Study area (km²):', waGeometry.area(1).divide(1e6));

// Generate all 4 environmental layers
var ndviLayer = computeNDVI(ee.Feature(waGeometry), startDate, endDate);
var no2Layer = computeNO2(ee.Feature(waGeometry), startDate, endDate);
var lstLayer = computeLST(ee.Feature(waGeometry), startDate, endDate);
var pm25Layer = computePM25(ee.Feature(waGeometry), startDate, endDate);

// Visualize environmental layers
Map.addLayer(ndviLayer, {
  min: 0.2, max: 0.8, 
  palette: ['white', 'yellow', 'green']
}, 'NDVI (Green Space)', false);

Map.addLayer(no2Layer, {
  min: 0, max: 0.0002,
  palette: ['blue', 'yellow', 'red']
}, 'NO2 (Air Pollution)');

Map.addLayer(lstLayer, {
  min: 10, max: 30,
  palette: ['blue', 'white', 'red']
}, 'Temperature (°C)', false);

Map.addLayer(pm25Layer, {
  min: 0, max: 30,
  palette: ['green', 'yellow', 'orange', 'red', 'purple']
}, 'PM2.5 (μg/m³)', false);

print('Environmental layers computed successfully!');
print('Toggle layers on/off in the Layers panel on the map.');

// Export NDVI
Export.image.toDrive({
  image: ndviLayer,
  description: 'WA_NDVI_2023_MayJune',
  scale: 30,
  region: waGeometry,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// Export NO2
Export.image.toDrive({
  image: no2Layer,
  description: 'WA_NO2_2023_MayJune',
  scale: 1000, // Coarser resolution for NO2
  region: waGeometry,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// Export Temperature
Export.image.toDrive({
  image: lstLayer,
  description: 'WA_LST_2023_MayJune',
  scale: 1000, // 1km MODIS resolution
  region: waGeometry,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// Export PM2.5
Export.image.toDrive({
  image: pm25Layer,
  description: 'WA_PM25_2023_MayJune',
  scale: 250, // GEOS-CF is ~25km, but save at finer scale
  region: waGeometry,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// Function to monitor computation progress
function getComputationSize(geometry) {
  var area = geometry.area(1).divide(1e6); // km²
  print('Processing area (km²):', area);
  return area;
}

// Legend helper with show/hide capability
function createColorLegend(title, palette, min, max, position) {
  var legendTitle = ui.Label(title, {fontWeight: 'bold', fontSize: '14px'});
  var colorBar = ui.Thumbnail({
    image: ee.Image.pixelLonLat().select(0)
      .multiply((max - min) / 100.0).add(min)
      .visualize({min: min, max: max, palette: palette}),
    params: {bbox: [0, 0, 100, 10], dimensions: '100x10'},
    style: {stretch: 'horizontal', margin: '0px 8px', maxHeight: '20px'}
  });
  var legendLabels = ui.Panel([
    ui.Label(min.toFixed(2)), 
    ui.Label(((min + max) / 2).toFixed(2)), 
    ui.Label(max.toFixed(2))
  ], ui.Panel.Layout.flow('horizontal'));
  var legendPanel = ui.Panel({
    widgets: [legendTitle, colorBar, legendLabels],
    style: {position: position || 'bottom-left', padding: '8px 15px'}
  });
  
  return legendPanel;
}

// Create all legends
var ndviLegend = createColorLegend('NDVI (Green Space)', ['white', 'yellow', 'green'], 0.2, 0.8, 'bottom-left');
var no2Legend = createColorLegend('NO2 (mol/m²)', ['blue', 'yellow', 'red'], 0, 0.0002, 'bottom-left');
var lstLegend = createColorLegend('Temperature (°C)', ['blue', 'white', 'red'], 10, 30, 'bottom-left');
var pm25Legend = createColorLegend('PM2.5 (μg/m³)', ['green', 'yellow', 'orange', 'red', 'purple'], 0, 30, 'bottom-left');

// Create control panel for legend toggles
var controlPanel = ui.Panel({
  widgets: [
    ui.Label('📊 LEGEND CONTROLS', {fontWeight: 'bold', fontSize: '16px', margin: '0 0 10px 0'}),
    ui.Button({
      label: 'Show NDVI Legend',
      onClick: function() {
        Map.remove(no2Legend);
        Map.remove(lstLegend);
        Map.remove(pm25Legend);
        Map.add(ndviLegend);
      }
    }),
    ui.Button({
      label: 'Show NO2 Legend',
      onClick: function() {
        Map.remove(ndviLegend);
        Map.remove(lstLegend);
        Map.remove(pm25Legend);
        Map.add(no2Legend);
      }
    }),
    ui.Button({
      label: 'Show Temperature Legend',
      onClick: function() {
        Map.remove(ndviLegend);
        Map.remove(no2Legend);
        Map.remove(pm25Legend);
        Map.add(lstLegend);
      }
    }),
    ui.Button({
      label: 'Show PM2.5 Legend',
      onClick: function() {
        Map.remove(ndviLegend);
        Map.remove(no2Legend);
        Map.remove(lstLegend);
        Map.add(pm25Legend);
      }
    }),
    ui.Button({
      label: 'Hide All Legends',
      onClick: function() {
        Map.remove(ndviLegend);
        Map.remove(no2Legend);
        Map.remove(lstLegend);
        Map.remove(pm25Legend);
      },
      style: {color: 'red'}
    })
  ],
  style: {position: 'top-right', padding: '8px', backgroundColor: 'white'}
});

// Add control panel to map
Map.add(controlPanel);

// Calculate statistics for each environmental layer across WA
var ndviStats = ndviLayer.reduceRegion({
  reducer: ee.Reducer.mean().combine({
    reducer2: ee.Reducer.minMax(),
    sharedInputs: true
  }),
  geometry: waGeometry,
  scale: 1000,
  maxPixels: 1e13
});

var no2Stats = no2Layer.reduceRegion({
  reducer: ee.Reducer.mean().combine({
    reducer2: ee.Reducer.minMax(),
    sharedInputs: true
  }),
  geometry: waGeometry,
  scale: 1000,
  maxPixels: 1e13
});

var lstStats = lstLayer.reduceRegion({
  reducer: ee.Reducer.mean().combine({
    reducer2: ee.Reducer.minMax(),
    sharedInputs: true
  }),
  geometry: waGeometry,
  scale: 1000,
  maxPixels: 1e13
});

var pm25Stats = pm25Layer.reduceRegion({
  reducer: ee.Reducer.mean().combine({
    reducer2: ee.Reducer.minMax(),
    sharedInputs: true
  }),
  geometry: waGeometry,
  scale: 1000,
  maxPixels: 1e13
});

print('ENVIRONMENTAL LAYER STATISTICS');
print('NDVI (Green Space) Stats:', ndviStats);
print('NO2 (Air Pollution) Stats:', no2Stats);
print('Temperature (°C) Stats:', lstStats);
print('PM2.5 Proxy (AOD) Stats:', pm25Stats);

/* 
 == TRIP LOGIC SECTION ==
- note: each trip is currently labeled as 'trip_id'
- to run just one, uncomment the desired trip and comment out others
- otherwise new trip setup must be implemented by simlply adding 
- another 'trip' such as 'trip_id2' fully to the implementation.
*/

/*
var trip = ee.Feature(
  ee.Geometry.LineString([
    [-117.1540, 46.7324],   // Pullman
    [-117.4260, 47.6588]    // Spokane
  ]),
  {'trip_id': 'Pullman_to_Spokane_2_vert'}
);
*/

/*
var trip = ee.Feature(
  ee.Geometry.LineString([
    [-117.1863, 46.7483],   // PHS
    [-117.1628, 46.7313]    // Terrell Library
  ]),
  {'trip_id': 'Pullman_High_to_WSU'}
);
*/

//US 195 approx
var trip = ee.Feature(
  ee.Geometry.LineString([
  [-117.1818, 46.7313],  // Pullman
  [-117.2952, 46.7815],  // Halfway Pullman-Colfax
  [-117.3644, 46.8802],  // Colfax
  [-117.3810, 47.0630],  // Steptoe
  [-117.3660, 47.1300],  // Thornton (approx)
  [-117.3802, 47.4289],  // Spangle
  [-117.4260, 47.6588]   // Spokane
]),
  { 'trip_id': 'Pullman_to_Spokane_7_vert' }
);
var geom = trip.geometry();
var distance_m = geom.length({'maxError': 1});
print('Distance (km):', distance_m.divide(1000));


// Add trip to map
Map.addLayer(trip, {color: 'black', width: 3}, 'Sample Trip');

function createPointsAlongLine(lineGeometry, intervalMeters) {
  var length = lineGeometry.length();
  var distances = ee.List.sequence(0, length, intervalMeters);

  var points = distances.map(function(dist) {
    var pt = lineGeometry.cutLines([dist]).geometries().get(0);
    return ee.Feature(ee.Geometry(pt));
  });

  return ee.FeatureCollection(points);
}

function calculateDetailedExposure(tripFeature, envLayer, samplingInterval) {
  var points = createPointsAlongLine(tripFeature.geometry(), samplingInterval);

  var sampledPoints = points.map(function(point) {
    var value = envLayer.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: point.geometry().buffer(25),
      scale: 250,
      maxPixels: 1e9
    });

    return point.set('exposure', value.values().get(0));
  });

  return tripFeature.set({
    'mean_exposure': sampledPoints.aggregate_mean('exposure'),
    'max_exposure': sampledPoints.aggregate_max('exposure'),
    'min_exposure': sampledPoints.aggregate_min('exposure'),
    'sampling_points_count': sampledPoints.size(),
    'sampling_points_fc': sampledPoints
  });
}

function calculateMultiExposure(tripFeature, layersDict, samplingInterval) {
  var tripWithExposures = tripFeature;

  var layerNames = ee.List(Object.keys(layersDict));

  // Loop over each layer
  layerNames.getInfo().forEach(function(name) {
    var safeName = name.replace(/\./g, ''); 
    var layer = layersDict[name];
    var exposure = calculateDetailedExposure(tripWithExposures, layer, samplingInterval);

    tripWithExposures = tripWithExposures
      .set('mean_' + safeName, exposure.get('mean_exposure'))
      .set('max_' + safeName, exposure.get('max_exposure'))
      .set('min_' + safeName, exposure.get('min_exposure'))
      .set('sampling_points_count_' + safeName, exposure.get('sampling_points_count'))
      .set('sampling_points_fc_' + safeName, exposure.get('sampling_points_fc'));
  });

  return tripWithExposures;
}

// == MultiExposure Calculation ==

// Dictionary of environmental layers
var layers = {
  'NDVI': ndviLayer,
  'NO2': no2Layer,
  'Temperature': lstLayer,
  'PM25': pm25Layer 
};

// Sampling interval in meters
var samplingInterval = 1000;

// Compute exposure
var multiExposure = calculateMultiExposure(trip, layers, samplingInterval);

print('Trip Exposure Summary (All Layers):', multiExposure);

// List of layer names 
var layerNames = ee.List(['NDVI', 'NO2', 'Temperature', 'PM25']);

// Build a FeatureCollection of {layer, stat, value} using properties on multiExposure
var statsFc = ee.FeatureCollection(
  layerNames.map(function(name) {
    name = ee.String(name);
    var meanVal = multiExposure.get(ee.String('mean_').cat(name));
    var minVal = multiExposure.get(ee.String('min_').cat(name));
    var maxVal = multiExposure.get(ee.String('max_').cat(name));
    // Create three features per layer: Min, Mean, Max
    return ee.List([
      ee.Feature(null, {layer: name, stat: 'Min',   value: minVal}),
      ee.Feature(null, {layer: name, stat: 'Mean',  value: meanVal}),
      ee.Feature(null, {layer: name, stat: 'Max',   value: maxVal})
    ]);
  }).flatten()
);

// Print the table in console
print('Trip stats table (layer, stat, value):', statsFc);

//small charts
layerNames.getInfo().forEach(function(name) {
  var layerFilter = statsFc.filter(ee.Filter.eq('layer', name));
  var smallChart = ui.Chart.feature.byFeature(layerFilter, 'stat', 'value')
    .setChartType('ColumnChart')
    .setOptions({
      title: name + ' — Min / Mean / Max',
      legend: {position: 'none'},
      hAxis: {title: ''},
      vAxis: {title: 'Value'},
    });
  print(smallChart);
});

