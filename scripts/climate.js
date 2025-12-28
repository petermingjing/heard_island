// Heard Island Climate Metrics to Match Albedo Variations
// Google Earth Engine JavaScript
// =============================================================================

// 1) STUDY AREA AND TEMPORAL SETUP
// =============================================================================
var heardIsland = ee.FeatureCollection('projects/ee-glacier/assets/Individual_glaciers/heard_island');
var glacierGeometry = heardIsland.geometry().dissolve();

var startDate = '2000-01-01';
var endDate   = '2025-01-01';

print('=== HIG CLIMATE METRICS (Daily, ERA5-Land) ===');
print('Glacier area (km²):', glacierGeometry.area().divide(1e6));
print('Time period:', startDate, 'to', endDate);

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

// Helper: return first-day-of-month date string
function monthKey(img) {
  var d = ee.Date(img.get('system:time_start'));
  return d.format('YYYY-MM-01');
}

// Helper: number of days in month for an image
function daysInMonth(img) {
  var d0 = ee.Date(img.get('system:time_start'));
  var d1 = d0.advance(1, 'month');
  return d1.difference(d0, 'day');
}

// 2) ERA5-Land datasets (HOURLY→daily): Temperature, Precipitation, Snowfall, SSRD
// =============================================================================
// Bands (typical): temperature_2m [K], total_precipitation [m], snowfall [m]
var era5Monthly = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY')
  .filterDate(startDate, endDate)
  .filterBounds(glacierGeometry);

print('ERA5-Land monthly images (context):', era5Monthly.size());

// Also set up ERA5-Land hourly for daily aggregation
var era5Hourly = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY')
  .filterDate(startDate, endDate)
  .filterBounds(glacierGeometry);
print('ERA5-Land hourly images:', era5Hourly.size());

// 3) DATE RANGES AND LIST BUILDERS
// =============================================================================
// Months list is optional; we mainly use the daily list below
var start = ee.Date(startDate);
var end = ee.Date(endDate);
var nMonths = end.difference(start, 'month').floor();
var monthsList = ee.List.sequence(0, nMonths.subtract(1)).map(function(i) {
  return start.advance(ee.Number(i), 'month');
});

// 4) BUILD DAILY SERIES (ERA5 only) WITHOUT HEAVY JOINS
// =============================================================================
var nDays = end.difference(start, 'day').floor();
var daysList = ee.List.sequence(0, nDays.subtract(1)).map(function(i) {
  return start.advance(ee.Number(i), 'day');
});

var dailySeries = ee.FeatureCollection(daysList.map(function(d) {
  d = ee.Date(d);
  var dNext = d.advance(1, 'day');

  var dayCol = era5Hourly.filterDate(d, dNext);
  var hasData = dayCol.size().gt(0);

  return ee.Feature(ee.Algorithms.If(hasData,
    (function() {
      var t2mDailyC = dayCol.select('temperature_2m').mean().subtract(273.15);
      var tpDaily_m = dayCol.select('total_precipitation').sum();
      var sfDaily_m = dayCol.select('snowfall').sum();
      var ssrdDaily_J = dayCol.select('surface_solar_radiation_downwards').sum();

      // Use a single-band reference projection to avoid mixed-projection errors
      var refProj = dayCol.select('temperature_2m').first().projection();
      var scale = refProj.nominalScale();
      var pixelArea = ee.Image.pixelArea().reproject(refProj);

      // Area-weighted mean helper
      function areaWeightedMean(img, band) {
        var imgRP = img.reproject(refProj);
        var num = imgRP.multiply(pixelArea).reduceRegion({
          reducer: ee.Reducer.sum(), geometry: glacierGeometry, scale: scale,
          maxPixels: 1e9, bestEffort: true, tileScale: 4
        }).get(band);
        var den = pixelArea.updateMask(imgRP.mask()).reduceRegion({
          reducer: ee.Reducer.sum(), geometry: glacierGeometry, scale: scale,
          maxPixels: 1e9, bestEffort: true, tileScale: 4
        }).get('area');
        return ee.Number(num).divide(ee.Number(den));
      }

      var tC = areaWeightedMean(t2mDailyC, 'temperature_2m');
      var ppt_m = areaWeightedMean(tpDaily_m, 'total_precipitation');
      var snow_m = areaWeightedMean(sfDaily_m, 'snowfall');

      // Area-sum volumes (m^3)
      var ppt_vol_m3 = tpDaily_m.reproject(refProj).multiply(pixelArea).reduceRegion({
        reducer: ee.Reducer.sum(), geometry: glacierGeometry, scale: scale,
        maxPixels: 1e9, bestEffort: true, tileScale: 4
      }).get('total_precipitation');
      var snow_vol_m3 = sfDaily_m.reproject(refProj).multiply(pixelArea).reduceRegion({
        reducer: ee.Reducer.sum(), geometry: glacierGeometry, scale: scale,
        maxPixels: 1e9, bestEffort: true, tileScale: 4
      }).get('snowfall');

      // SSRD area-weighted mean (J/m^2 per day) -> W/m^2
      var ssrd_J = areaWeightedMean(ssrdDaily_J, 'surface_solar_radiation_downwards');
      var ssrd_wm2 = ee.Number(ssrd_J).divide(86400);

      var ppt_mm = ee.Number(ppt_m).multiply(1000);
      var snow_mm = ee.Number(snow_m).multiply(1000);
      var pdd = ee.Number(0).max(ee.Number(tC));

      return ee.Feature(null, {
        date: d.format('YYYY-MM-dd'),
        timestamp: d.millis(),
        source: 'ERA5-Land Daily',
        t2m_C: tC,
        pdd_degCday: pdd,
        precip_mm: ppt_mm,
        snowfall_mm_we: snow_mm,
        precip_volume_m3: ppt_vol_m3,
        snowfall_volume_m3: snow_vol_m3,
        ssrd_wm2: ssrd_wm2
      });
    })(),
    ee.Feature(null, {
      date: d.format('YYYY-MM-dd'),
      timestamp: d.millis(),
      source: 'ERA5-Land Daily',
      t2m_C: null,
      pdd_degCday: null,
      precip_mm: null,
      snowfall_mm_we: null,
      precip_volume_m3: null,
      snowfall_volume_m3: null,
      ssrd_wm2: null
    })
  ));
}));

// 5) SUMMARY (console-safe)
// =============================================================================
// Avoid printing the full daily FeatureCollection; just print expected row count
print('Daily ERA5 climate rows expected:', nDays);

// 6) OPTIONAL VISUALIZATION
// =============================================================================
Map.centerObject(glacierGeometry, 9);
Map.addLayer(glacierGeometry, {color: 'orange'}, 'Heard Island Glaciers');

// Optional: skip heavy daily charts to avoid console limits

// 7) EXPORTS
// =============================================================================
Export.table.toDrive({
  collection: ee.FeatureCollection(dailySeries),
  description: 'HIG_ERA5_Climate_Daily',
  folder: 'HIG_climate_exports',
  fileNamePrefix: 'heard_island_era5_land_climate_daily',
  fileFormat: 'CSV'
});

print('=== EXPORT QUEUED ===');
print('ERA5 climate daily CSV -> Drive/HIG_climate_exports');