// Heard Island SRTM DEM Export
// Google Earth Engine Code Editor
// Simple script to export SRTM elevation data for Heard Island
// =============================================================================

// 1. STUDY AREA SETUP
// =============================================================================
// Define bounding box for entire Heard Island
// Island center: 53°06'S, 73°31'E
// Island covers ~368 km², approximately 40 km long
// Bounding box: roughly 53.0°S to 53.3°S, 73.0°E to 73.7°E
var islandBounds = ee.Geometry.Rectangle([
  73.0,  // West longitude
  -53.3, // South latitude
  73.7,  // East longitude
  -53.0  // North latitude
]);

// Optional: Also load glacier geometry for reference
var heardIslandGlaciers = ee.FeatureCollection("projects/ee-glacier/assets/Individual_glaciers/heard_island");
var glacierGeometry = heardIslandGlaciers.geometry().dissolve();

print('=== HEARD ISLAND SRTM DEM EXPORT ===');
print('Island bounding box area (km²):', islandBounds.area().divide(1e6));
print('Glacier area (km²):', glacierGeometry.area().divide(1e6));
print('');

// 2. LOAD SRTM DEM
// =============================================================================
// USGS/SRTMGL1_003: 30m resolution, collected February 2000
var srtmDEM = ee.Image("USGS/SRTMGL1_003").clip(islandBounds);

// Get elevation statistics for entire island
var elevationStats = srtmDEM.reduceRegion({
  reducer: ee.Reducer.minMax().combine({
    reducer2: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    sharedInputs: true
  }),
  geometry: islandBounds,
  scale: 30,
  maxPixels: 1e9,
  bestEffort: true
});

print('SRTM Elevation Statistics (entire island):');
print('  Min elevation (m):', elevationStats.get('elevation_min'));
print('  Max elevation (m):', elevationStats.get('elevation_max'));
print('  Mean elevation (m):', elevationStats.get('elevation_mean'));
print('  StdDev elevation (m):', elevationStats.get('elevation_stdDev'));
print('');

// 3. VISUALIZATION
// =============================================================================
Map.centerObject(islandBounds, 9);
Map.addLayer(islandBounds, {color: 'blue', opacity: 0.2}, 'Heard Island Bounds');
Map.addLayer(glacierGeometry, {color: 'orange'}, 'Heard Island Glaciers (reference)');
Map.addLayer(srtmDEM, {
  min: 0,
  max: 2745,
  palette: ['blue', 'green', 'brown', 'white']
}, 'SRTM Elevation');

// 4. EXPORT SRTM DEM
// =============================================================================
Export.image.toDrive({
  image: srtmDEM,
  description: 'HIG_SRTM_DEM',
  folder: 'HIG_SRTM_exports',
  fileNamePrefix: 'heard_island_srtm_dem',
  region: islandBounds,
  scale: 30,                    // 30m resolution
  crs: 'EPSG:4326',            // WGS84 geographic coordinate system
  maxPixels: 1e9,              // Maximum pixels allowed
  formatOptions: {
    cloudOptimized: true       // Create cloud-optimized GeoTIFF
  }
});

print('=== EXPORT QUEUED ===');
print('SRTM DEM will be exported as GeoTIFF to:');
print('  Google Drive → HIG_SRTM_exports → heard_island_srtm_dem');
print('');
print('Check the Tasks tab (right panel) to start the export.');
print('The export will create a 30m resolution DEM covering Heard Island.');
print('');
print('Note: SRTM data is from February 2000 (Shuttle Radar Topography Mission)');

