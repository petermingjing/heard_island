// Heard Island Mass Balance Data Extraction
// Google Earth Engine Code Editor
// Simple script to extract mass balance data from multiple sources
// =============================================================================

// 1. STUDY AREA SETUP
// =============================================================================
var heardIsland = ee.FeatureCollection("projects/ee-glacier/assets/Individual_glaciers/heard_island");
var glacierGeometry = heardIsland.geometry().dissolve();

// Temporal parameters
var startDate = '2000-01-01';
var endDate = '2025-01-01';

print('=== HEARD ISLAND MASS BALANCE DATA EXTRACTION ===');
print('Glacier area (km²):', glacierGeometry.area().divide(1e6));
print('Time period:', startDate, 'to', endDate);
print('');

// Console-safe printing helpers
var DEBUG = false; // set true to see up to 20 sample rows
function safePrintFC(label, fc) {
  var fcf = ee.FeatureCollection(fc);
  print(label + ' size:', ee.Number(fcf.size()));
  if (DEBUG) {
    print(label + ' sample (up to 20):', fcf.limit(20));
  } else {
    print(label + ' sample (1):', fcf.limit(1));
  }
}

// 2. DATA SOURCES
// =============================================================================

// GRACE Mass Balance Data
var graceCollection = ee.ImageCollection("NASA/GRACE/MASS_GRIDS_V04/MASCON_CRI")
  .filterDate(startDate, endDate)
  .filterBounds(glacierGeometry);

// SRTM DEM for geodetic analysis
var srtmDEM = ee.Image("USGS/SRTMGL1_003").clip(glacierGeometry);

// NASADEM (refined SRTM) as alternative 2000 baseline
var nasademDEM = ee.Image("NASA/NASADEM_HGT/001")
  .select('elevation')
  .clip(glacierGeometry);

// Additional DEMs for post-2000 comparison
var alosDEM = ee.ImageCollection("JAXA/ALOS/AW3D30/V3_2")
  .select('DSM')
  .mosaic()
  .clip(glacierGeometry);

var copernicusDEM = ee.ImageCollection("COPERNICUS/DEM/GLO30")
  .select('DEM')
  .mosaic()
  .clip(glacierGeometry);

print('Available data sources:');
print('- GRACE observations:', graceCollection.size());
print('');

// 3. GRACE MASS BALANCE EXTRACTION
// =============================================================================
function extractGRACE(image) {
  var lwe = image.select('lwe_thickness'); // cm water equivalent
  var uncertainty = image.select('uncertainty'); // cm (std)

  // Build fractional-overlap weights on the native GRACE grid
  var graceProj = image.projection();
  var hiMask = ee.Image.constant(1)
    .clip(glacierGeometry)
    .rename('mask')
    .reproject({crs: 'EPSG:3857', scale: 2000});
  var frac = hiMask
    .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 65535})
    .reproject(graceProj)
    .unmask(0);
  var pixelArea = ee.Image.pixelArea().reproject(graceProj);
  var wArea = frac.multiply(pixelArea); // m^2 overlap per mascon

  // Weighted mean: sum(lwe * wArea) / sum(wArea)
  var num = lwe.multiply(wArea);
  var sumNum = num.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: glacierGeometry,
    scale: graceProj.nominalScale(),
    maxPixels: 1e9,
    bestEffort: true,
    tileScale: 4
  }).get('lwe_thickness');
  var sumDen = wArea.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: glacierGeometry,
    scale: graceProj.nominalScale(),
    maxPixels: 1e9,
    bestEffort: true,
    tileScale: 4
  }).get('mask');

  var weightedMean = ee.Algorithms.If(
    ee.Algorithms.IsEqual(sumDen, null),
    null,
    ee.Number(sumNum).divide(ee.Number(sumDen))
  );

  // Weighted uncertainty (independent mascons):
  // sigma_w = sqrt( sum( (sigma_i^2 * w_i^2) ) ) with w_i = wArea_i / sumDen
  var wAreaSq = wArea.multiply(wArea);
  var uncSq = uncertainty.multiply(uncertainty);
  var uncTerm = uncSq.multiply(wAreaSq);
  var sumUncTerm = uncTerm.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: glacierGeometry,
    scale: graceProj.nominalScale(),
    maxPixels: 1e9,
    bestEffort: true,
    tileScale: 4
  }).get('uncertainty');
  var weightedSigma = ee.Algorithms.If(
    ee.Algorithms.IsEqual(sumDen, null),
    null,
    ee.Number(sumUncTerm).divide(ee.Number(sumDen).pow(2)).sqrt()
  );

  return ee.Feature(null, {
    'date': image.date().format('YYYY-MM-dd'),
    'timestamp': image.date().millis(),
    'source': 'GRACE',
    'mass_balance_cm': weightedMean,
    'uncertainty_cm': weightedSigma
  });
}

var graceTimeSeries = graceCollection.map(extractGRACE);

// 4. GEODETIC MASS BALANCE (Simple Elevation Extraction)
// =============================================================================
function extractElevation(image, source) {
  var elevation = image.select(0); // First band
  
  var stats = elevation.reduceRegion({
    reducer: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.stdDev(), 
      sharedInputs: true
    }),
    geometry: glacierGeometry,
    scale: 30,
    maxPixels: 1e9,
    bestEffort: true,
    tileScale: 4
  });
  
  return ee.Feature(null, {
    'date': image.date().format('YYYY-MM-dd'),
    'timestamp': image.date().millis(),
    'source': source,
    'elevation_mean': ee.Algorithms.If(
      stats.contains('elevation_mean'),
      stats.get('elevation_mean'),
      null
    ),
    'elevation_stddev': ee.Algorithms.If(
      stats.contains('elevation_stdDev'), 
      stats.get('elevation_stdDev'), 
      null
    )
  });
}

// Geodetic mass balance relative to SRTM baseline
function computeGeodeticMassBalance(laterDEM, laterDateString, source, baselineDEM, baselineName) {
  var dh = laterDEM.subtract(baselineDEM).rename('dh_m');
  
  var stats = dh.reduceRegion({
    reducer: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    geometry: glacierGeometry,
    scale: 30,
    maxPixels: 1e9,
    bestEffort: true,
    tileScale: 4
  });
  
  var rhoIce = 900;   // kg/m^3
  var rhoWater = 1000; // kg/m^3
  var dhMean = ee.Algorithms.If(stats.contains('dh_m_mean'), stats.get('dh_m_mean'), null);
  var dhStd = ee.Algorithms.If(stats.contains('dh_m_stdDev'), stats.get('dh_m_stdDev'), null);
  var weMeters = ee.Algorithms.If(dhMean, ee.Number(dhMean).multiply(rhoIce).divide(rhoWater), null);
  var weCm = ee.Algorithms.If(weMeters, ee.Number(weMeters).multiply(100), null);
  var weCmStd = ee.Algorithms.If(dhStd, ee.Number(dhStd).multiply(rhoIce).divide(rhoWater).multiply(100), null);
  
  var dateMs = ee.Date(laterDateString).millis();
  return ee.Feature(null, {
    'date': ee.Date(laterDateString).format('YYYY-MM-dd'),
    'timestamp': dateMs,
    'source': 'Geodetic_' + source,
    'baseline': baselineName,
    'dh_m': dhMean,
    'dh_std_m': dhStd,
    'mass_balance_cm_we': weCm,
    'mass_balance_cm_we_std': weCmStd,
    'density_assumed_kg_m3': rhoIce
  }).set('system:time_start', dateMs);
}

// SRTM baseline elevation
var srtmFeature = extractElevation(srtmDEM.set('system:time_start', ee.Date('2000-02-11').millis()), 'SRTM');

// NASADEM baseline elevation
var nasademFeature = extractElevation(nasademDEM.set('system:time_start', ee.Date('2000-02-11').millis()), 'NASADEM');

// Geodetic mass balance features relative to SRTM and NASADEM (approximate epochs)
var geodeticAlosSRTM = computeGeodeticMassBalance(
  alosDEM.set('system:time_start', ee.Date('2010-07-01').millis()),
  '2010-07-01',
  'ALOS_AW3D30',
  srtmDEM,
  'USGS/SRTMGL1_003'
);

var geodeticAlosNASADEM = computeGeodeticMassBalance(
  alosDEM.set('system:time_start', ee.Date('2010-07-01').millis()),
  '2010-07-01',
  'ALOS_AW3D30',
  nasademDEM,
  'NASA/NASADEM_HGT/001'
);

var geodeticCopernicusSRTM = computeGeodeticMassBalance(
  copernicusDEM.set('system:time_start', ee.Date('2020-07-01').millis()),
  '2020-07-01',
  'Copernicus_GLO30',
  srtmDEM,
  'USGS/SRTMGL1_003'
);

var geodeticCopernicusNASADEM = computeGeodeticMassBalance(
  copernicusDEM.set('system:time_start', ee.Date('2020-07-01').millis()),
  '2020-07-01',
  'Copernicus_GLO30',
  nasademDEM,
  'NASA/NASADEM_HGT/001'
);

var geodeticSeries = ee.FeatureCollection([
  geodeticAlosSRTM,
  geodeticAlosNASADEM,
  geodeticCopernicusSRTM,
  geodeticCopernicusNASADEM
]);

// 5. COMBINE ALL DATA
// =============================================================================
var allMassBalanceData = graceTimeSeries
  .merge(ee.FeatureCollection([srtmFeature, nasademFeature]))
  .merge(geodeticSeries);

// 6. SUMMARY STATISTICS
// =============================================================================
print('=== DATA EXTRACTION SUMMARY ===');
print('Total GRACE records:', graceTimeSeries.size());
print('Total baseline features (SRTM + NASADEM):', ee.Number(2));
print('Total Geodetic records:', geodeticSeries.size());
print('Total combined records:', allMassBalanceData.size());
print('');

// Sample of GRACE data (console-safe)
safePrintFC('GRACE time series', graceTimeSeries);

// 7. SIMPLE VISUALIZATION
// =============================================================================
Map.centerObject(glacierGeometry, 10);
Map.addLayer(glacierGeometry, {color: 'red'}, 'Heard Island Glaciers');
Map.addLayer(srtmDEM, {min: 0, max: 2745, palette: ['blue', 'green', 'brown', 'white']}, 'SRTM Elevation');
Map.addLayer(nasademDEM, {min: 0, max: 2745, palette: ['navy', 'green', 'brown', 'white']}, 'NASADEM Elevation');

// GRACE mascon pixels intersecting Heard Island
var graceSampleImage = graceCollection.first();
var graceProj = graceSampleImage.projection();
print('GRACE projection:', graceProj);

// Rasterize glacier at high resolution, then aggregate to GRACE grid to include any intersecting pixels
var highResMask = ee.Image.constant(1)
  .clip(glacierGeometry)
  .rename('mask')
  .reproject({crs: 'EPSG:3857', scale: 2000});

var glacierMaskOnGraceGrid = highResMask
  .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 65535})
  .reproject(graceProj)
  .gt(0)
  .selfMask();

// Vectorize those pixels for display
var masconPixels = glacierMaskOnGraceGrid.reduceToVectors({
  geometry: glacierGeometry.buffer(50000),
  scale: graceProj.nominalScale(),
  geometryType: 'polygon',
  labelProperty: 'mascon',
  maxPixels: 1e9
});

print('Mascon pixels intersecting glacier:', masconPixels.size());
Map.addLayer(ee.FeatureCollection(masconPixels), {color: '#FFFF00'}, 'GRACE Mascons Intersecting Glacier');

// Optional: show boundaries of each GRACE mascon cell around the glacier (200 km buffer)
var masconGridAOI = glacierGeometry.buffer(200000);
var px = ee.Image.pixelCoordinates(graceProj);
// Create a unique ID per pixel using (x, y)
var pxId = px.select('x').multiply(1000000).add(px.select('y')).toInt64();
var masconGridVals = pxId.clip(masconGridAOI);
var masconGridVectors = masconGridVals.reduceToVectors({
  geometry: masconGridAOI,
  scale: graceProj.nominalScale(),
  geometryType: 'polygon',
  labelProperty: 'px_id',
  maxPixels: 1e9
});
print('Mascon grid cells (buffered AOI):', masconGridVectors.size());
var masconGridOutline = ee.FeatureCollection(masconGridVectors).style({
  color: '#00FFFF',
  fillColor: '00000000',
  width: 1
});
Map.addLayer(masconGridOutline, {}, 'GRACE Mascon Grid (Outline)');

// Prepare intersecting mascons with pixel indices for export
var masconPixelsWithXY = ee.FeatureCollection(masconPixels).map(function(f) {
  var c = ee.Feature(f).centroid(1);
  var xy = ee.Image.pixelCoordinates(graceProj).sample(c.geometry(), graceProj.nominalScale()).first();
  return ee.Feature(f).set({
    'px_x': ee.Number(xy.get('x')),
    'px_y': ee.Number(xy.get('y'))
  });
});
print('Intersecting GRACE mascons (with indices):', masconPixelsWithXY.size());

// GRACE time series chart
var graceChart = ui.Chart.feature.byFeature(
  graceTimeSeries, 
  'date', 
  'mass_balance_cm'
)
.setChartType('LineChart')
.setOptions({
  title: 'GRACE Mass Balance Time Series (cm water equivalent)',
  hAxis: {title: 'Date'},
  vAxis: {title: 'Mass Balance (cm)'},
  lineWidth: 2,
  pointSize: 3
});
print(graceChart);

// Geodetic mass balance chart (force numeric x-axis) and ignore null values
var geodeticSeriesForChart = geodeticSeries
  .filter(ee.Filter.notNull(['mass_balance_cm_we', 'system:time_start']))
  .map(function(f) {
    return ee.Feature(f).set('timestamp_num', ee.Number(f.get('system:time_start')));
  });

print('Geodetic features available for chart (non-null):', geodeticSeriesForChart.size());

var geodeticChart = ui.Chart.feature.byFeature(
  geodeticSeriesForChart,
  'timestamp_num',
  'mass_balance_cm_we'
)
.setChartType('ColumnChart')
.setOptions({
  title: 'Geodetic Mass Balance (cm w.e., relative to SRTM/NASADEM)',
  hAxis: {title: 'Date (ms since epoch)'},
  vAxis: {title: 'Mass Balance (cm w.e.)'},
  legend: {position: 'none'}
});
print(geodeticChart);

// 8. DATA EXPORTS
// =============================================================================
// Export all mass balance data
Export.table.toDrive({
  collection: allMassBalanceData,
  description: 'HIG_all_mass_balance_data',
  folder: 'HIG_mass_balance_exports',
  fileNamePrefix: 'heard_island_mass_balance_all_sources',
  fileFormat: 'CSV'
});

// Export GRACE data only
Export.table.toDrive({
  collection: graceTimeSeries,
  description: 'HIG_GRACE_data',
  folder: 'HIG_mass_balance_exports',
  fileNamePrefix: 'heard_island_GRACE_mass_balance',
  fileFormat: 'CSV'
});

// Export baseline elevation features (SRTM + NASADEM)
var baselineElevationData = ee.FeatureCollection([srtmFeature, nasademFeature]);

Export.table.toDrive({
  collection: baselineElevationData,
  description: 'HIG_baseline_elevation_features',
  folder: 'HIG_mass_balance_exports',
  fileNamePrefix: 'heard_island_baseline_elevation_features',
  fileFormat: 'CSV'
});

// Export geodetic mass balance
Export.table.toDrive({
  collection: geodeticSeries,
  description: 'HIG_geodetic_mass_balance',
  folder: 'HIG_mass_balance_exports',
  fileNamePrefix: 'heard_island_geodetic_mass_balance',
  fileFormat: 'CSV'
});

// Export SRTM DEM
Export.image.toDrive({
  image: srtmDEM,
  description: 'HIG_SRTM_DEM',
  folder: 'HIG_mass_balance_exports',
  fileNamePrefix: 'heard_island_srtm_dem',
  region: glacierGeometry,
  scale: 30,
  crs: 'EPSG:4326',
  maxPixels: 1e9
});

// Export NASADEM DEM
Export.image.toDrive({
  image: nasademDEM,
  description: 'HIG_NASADEM_DEM',
  folder: 'HIG_mass_balance_exports',
  fileNamePrefix: 'heard_island_nasadem_dem',
  region: glacierGeometry,
  scale: 30,
  crs: 'EPSG:4326',
  maxPixels: 1e9
});

// Export intersecting GRACE mascons as Shapefile
Export.table.toDrive({
  collection: masconPixelsWithXY,
  description: 'HIG_GRACE_Mascons_Intersecting',
  folder: 'HIG_mass_balance_exports',
  fileNamePrefix: 'heard_island_grace_mascons_intersecting',
  fileFormat: 'SHP'
});

print('=== EXPORTS QUEUED ===');
print('Check the Tasks tab for:');
print('1. All mass balance data (combined CSV)');
print('2. GRACE data only (CSV)');
print('3. Baseline elevation features (CSV)');
print('4. SRTM DEM (GeoTIFF)');
print('5. NASADEM DEM (GeoTIFF)');
print('6. Geodetic mass balance (CSV)');
print('');
print('=== DATA USAGE NOTES ===');
print('- GRACE: Regional mass balance trends (25km resolution)');
print('- SRTM/NASADEM: Baseline elevation (2000, 30m resolution)');
print('- ALOS & Copernicus DEMs: Later epochs for geodetic differencing');
print('');
print('⚠️  GRACE resolution is too coarse for Heard Island (~311 km²)');
print('📊 Use GRACE for regional context; use DEM differencing for local geodetic mass balance');
