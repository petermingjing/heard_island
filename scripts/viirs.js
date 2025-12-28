// Heard Island Glacier VIIRS Albedo Analysis
// Google Earth Engine JavaScript
// =============================================================================
//
// Datasets:
// - NOAA/VIIRS/001/VNP43IA1 (BRDF/Albedo Model Parameters) - for albedo calculation
// - NOAA/VIIRS/001/VNP43IA2 (BRDF/Albedo Quality) - for quality assessment
// Quality values: 0=best (full inversion), 1=good, 2=fair, 3=magnitude inversion
//
// =============================================================================

// 1. STUDY AREA SETUP
// =============================================================================
var heardIsland = ee.FeatureCollection("projects/ee-glacier/assets/Individual_glaciers/heard_island");
var glacierGeometry = heardIsland.geometry();
// Console-safe printing helpers
var DEBUG = false; // set true to see up to 20 sample items
function safePrintIC(label, ic) {
  var icc = ee.ImageCollection(ic);
  print(label + ' size:', icc.size());
  var first = icc.first();
  print(label + ' first date:', ee.Algorithms.If(first,
    ee.Date(first.get('system:time_start')).format('YYYY-MM-dd'), 'NA'));
  if (DEBUG) {
    print(label + ' sample (1 image):', first);
  }
}

// Create interior mask for VIIRS (375m resolution)
var viirsProj = ee.Projection('EPSG:4326').atScale(375);

var interiorMask = ee.Image.constant(1)
  .clip(glacierGeometry)
  .reproject({crs: viirsProj})
  .unmask(0)
  .focal_min(1, 'circle')
  .rename('interior_mask');

var startDate = '2012-01-18';
var endDate = '2025-09-01';
var driveFolder = 'HIG_VIIRS_Albedo';

// =============================================================================
// Collections
// =============================================================================
var viirsBRDF = ee.ImageCollection('NOAA/VIIRS/001/VNP43IA1')
  .filterDate(startDate, endDate)
  .filterBounds(glacierGeometry);

var viirsQuality = ee.ImageCollection('NOAA/VIIRS/001/VNP43IA2')
  .filterDate(startDate, endDate)
  .filterBounds(glacierGeometry);

// =============================================================================
// 2. VNP43IA1 BRDF PROCESSING FUNCTION
// =============================================================================

function prepareVIIRSBRDF(image) {
  // VNP43IA1 provides BRDF/Albedo Model Parameters
  // Note: Parameters in GEE appear already scaled; avoid multiplying by 0.001 again
  
  // Per-band mandatory quality (0=best,1=good)
  var qI1 = image.select('BRDF_Albedo_Band_Mandatory_Quality_I1').lte(1);
  var qI2 = image.select('BRDF_Albedo_Band_Mandatory_Quality_I2').lte(1);
  var qI3 = image.select('BRDF_Albedo_Band_Mandatory_Quality_I3').lte(1);

  // BRDF parameters for I1 band (Blue)
  var fisoI1 = image.select('BRDF_Albedo_Parameters_fiso_I1').updateMask(qI1).rename('fiso_i1');
  var fvolI1 = image.select('BRDF_Albedo_Parameters_fvol_I1').updateMask(qI1).rename('fvol_i1');
  var fgeoI1 = image.select('BRDF_Albedo_Parameters_fgeo_I1').updateMask(qI1).rename('fgeo_i1');
  
  // BRDF parameters for I2 band (Green)
  var fisoI2 = image.select('BRDF_Albedo_Parameters_fiso_I2').updateMask(qI2).rename('fiso_i2');
  var fvolI2 = image.select('BRDF_Albedo_Parameters_fvol_I2').updateMask(qI2).rename('fvol_i2');
  var fgeoI2 = image.select('BRDF_Albedo_Parameters_fgeo_I2').updateMask(qI2).rename('fgeo_i2');
  
  // BRDF parameters for I3 band (Red)
  var fisoI3 = image.select('BRDF_Albedo_Parameters_fiso_I3').updateMask(qI3).rename('fiso_i3');
  var fvolI3 = image.select('BRDF_Albedo_Parameters_fvol_I3').updateMask(qI3).rename('fvol_i3');
  var fgeoI3 = image.select('BRDF_Albedo_Parameters_fgeo_I3').updateMask(qI3).rename('fgeo_i3');
  
  var out = fisoI1.addBands(fvolI1).addBands(fgeoI1)
    .addBands(fisoI2).addBands(fvolI2).addBands(fgeoI2)
    .addBands(fisoI3).addBands(fvolI3).addBands(fgeoI3)
    .updateMask(interiorMask);
    
  return out.copyProperties(image, ['system:time_start', 'system:time_end']);
}

// =============================================================================
// 3. VNP43IA2 QUALITY PROCESSING FUNCTION
// =============================================================================

function prepareVIIRSQuality(image) {
  // VNP43IA2 provides BRDF/Albedo Band Quality flags
  // Quality values: 0=best (full inversion), 1=good, 2=fair, 3=magnitude inversion
  
  var qualityI1 = image.select('BRDF_Albedo_Band_Quality_I1').rename('quality_i1');
  var qualityI2 = image.select('BRDF_Albedo_Band_Quality_I2').rename('quality_i2');
  var qualityI3 = image.select('BRDF_Albedo_Band_Quality_I3').rename('quality_i3');
  
  // Days of valid observation for each band
  var validDaysI1 = image.select('BRDF_Albedo_ValidObs_I1').rename('valid_days_i1');
  var validDaysI2 = image.select('BRDF_Albedo_ValidObs_I2').rename('valid_days_i2');
  var validDaysI3 = image.select('BRDF_Albedo_ValidObs_I3').rename('valid_days_i3');
  
  // Snow flag (0=snow-free, 1=snow present)
  var snowFlag = image.select('Snow_BRDF_Albedo').rename('snow_flag');
  
  // Land/Water type classification
  var landWaterType = image.select('BRDF_Albedo_LandWaterType').rename('land_water_type');
  
  var out = qualityI1.addBands(qualityI2)
    .addBands(qualityI3)
    .addBands(validDaysI1)
    .addBands(validDaysI2)
    .addBands(validDaysI3)
    .addBands(snowFlag)
    .addBands(landWaterType)
    .updateMask(interiorMask);
    
  return out.copyProperties(image, ['system:time_start', 'system:time_end']);
}

// =============================================================================
// 4. VIIRS COLLECTIONS SETUP
// =============================================================================

// VNP43IA1 - BRDF/Albedo Model Parameters
var viirsBRDF = ee.ImageCollection('NOAA/VIIRS/001/VNP43IA1')
  .filterDate(startDate, endDate)
  .filterBounds(glacierGeometry);

var viirsBRDFProcessed = viirsBRDF.map(prepareVIIRSBRDF);

// VNP43IA2 - BRDF/Albedo Quality
var viirsQuality = ee.ImageCollection('NOAA/VIIRS/001/VNP43IA2')
  .filterDate(startDate, endDate)
  .filterBounds(glacierGeometry);

var viirsQualityProcessed = viirsQuality.map(prepareVIIRSQuality);

// =============================================================================
// 5. ALBEDO CALCULATION FUNCTION
// =============================================================================

function calculateAlbedoFromBRDF(image) {
  var img = ee.Image(image);
  // Calculate White-Sky Albedo (WSA) from BRDF parameters
  // Formula: WSA = fiso + 0.189184 * fvol - 1.377622 * fgeo
  
  // Albedo for I1 band (Blue)
  var albedoI1 = img.expression(
    'fiso + 0.189184 * fvol - 1.377622 * fgeo',
    {
      'fiso': img.select('fiso_i1'),
      'fvol': img.select('fvol_i1'),
      'fgeo': img.select('fgeo_i1')
    }
  ).rename('albedo_i1').clamp(0, 1);
  
  // Albedo for I2 band (Green)
  var albedoI2 = img.expression(
    'fiso + 0.189184 * fvol - 1.377622 * fgeo',
    {
      'fiso': img.select('fiso_i2'),
      'fvol': img.select('fvol_i2'),
      'fgeo': img.select('fgeo_i2')
    }
  ).rename('albedo_i2').clamp(0, 1);
  
  // Albedo for I3 band (Red)
  var albedoI3 = img.expression(
    'fiso + 0.189184 * fvol - 1.377622 * fgeo',
    {
      'fiso': img.select('fiso_i3'),
      'fvol': img.select('fvol_i3'),
      'fgeo': img.select('fgeo_i3')
    }
  ).rename('albedo_i3').clamp(0, 1);
  
  // Calculate broadband albedo (average of I1, I2, I3)
  var broadbandAlbedo = albedoI1.add(albedoI2).add(albedoI3).divide(3).rename('albedo_broadband');
  
  return albedoI1.addBands(albedoI2).addBands(albedoI3).addBands(broadbandAlbedo);
}

// Apply albedo calculation to BRDF data
var viirsAlbedoProcessed = viirsBRDFProcessed.map(calculateAlbedoFromBRDF);
// Preserve time properties on albedo images for charts/exports
viirsAlbedoProcessed = viirsBRDF.map(function(img) {
  var albedo = calculateAlbedoFromBRDF(prepareVIIRSBRDF(img));
  return albedo.copyProperties(img, ['system:time_start', 'system:time_end']);
});

// =============================================================================
// 6. VISUALIZATION
// =============================================================================
Map.centerObject(glacierGeometry, 9);

// Albedo visualization parameters
var albedoVis = {
  min: 0.0,
  max: 0.9,
  palette: ['#1a1a1a', '#c0c0c0', '#ffffff']
};

// Quality visualization parameters
var qualityVis = {
  min: 0,
  max: 3,
  palette: ['#00FF00', '#FFFF00', '#FF8000', '#FF0000'] // Green=best, Red=worst
};

var snowVis = {
  min: 0,
  max: 1,
  palette: ['#000080', '#FFFFFF'] // Blue=snow-free, White=snow
};

// Add albedo layers
Map.addLayer(viirsAlbedoProcessed.median().select('albedo_i1'), albedoVis, 'VIIRS Albedo I1');
Map.addLayer(viirsAlbedoProcessed.median().select('albedo_i2'), albedoVis, 'VIIRS Albedo I2');
Map.addLayer(viirsAlbedoProcessed.median().select('albedo_i3'), albedoVis, 'VIIRS Albedo I3');
Map.addLayer(viirsAlbedoProcessed.median().select('albedo_broadband'), albedoVis, 'VIIRS Broadband Albedo');

// Add quality layers
Map.addLayer(viirsQualityProcessed.median().select('quality_i1'), qualityVis, 'VIIRS Quality I1');
Map.addLayer(viirsQualityProcessed.median().select('snow_flag'), snowVis, 'VIIRS Snow Flag');

var chartAlbedoBroadband = ui.Chart.image.series({
  imageCollection: viirsAlbedoProcessed.select('albedo_broadband'),
  region: glacierGeometry,
  reducer: ee.Reducer.mean(),
  scale: 375
}).setOptions({
  title: 'VIIRS Broadband Albedo Time Series',
  hAxis: { title: 'Date' },
  vAxis: { title: 'Albedo', viewWindow: {min: 0, max: 1} },
  legend: { position: 'none' }
});
print('VIIRS Broadband Albedo Chart:', chartAlbedoBroadband);

var chartSnowFlag = ui.Chart.image.series({
  imageCollection: viirsQualityProcessed.select('snow_flag'),
  region: glacierGeometry,
  reducer: ee.Reducer.mean(),
  scale: 375
}).setOptions({
  title: 'VIIRS Snow Flag Time Series',
  hAxis: { title: 'Date' },
  vAxis: { title: 'Snow Flag', viewWindow: {min: 0, max: 1} },
  legend: { position: 'none' }
});
print('VIIRS Snow Flag Chart:', chartSnowFlag);

// =============================================================================
// 5. AREA-MEAN EXPORTS (CSV)
// =============================================================================
var areaMeanQuality = viirsQualityProcessed.map(function(img) {
  var meanDict = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: glacierGeometry,
    scale: 375,
    maxPixels: 1e13,
    tileScale: 4
  });
  
  return ee.Feature(null, {
    date: ee.Date(img.get('system:time_start')).format('YYYY-MM-dd'),
    quality_i1: meanDict.get('quality_i1'),
    quality_i2: meanDict.get('quality_i2'),
    quality_i3: meanDict.get('quality_i3'),
    valid_days_i1: meanDict.get('valid_days_i1'),
    valid_days_i2: meanDict.get('valid_days_i2'),
    valid_days_i3: meanDict.get('valid_days_i3'),
    snow_flag: meanDict.get('snow_flag'),
    land_water_type: meanDict.get('land_water_type'),
    product: 'VNP43IA2',
    source: 'VIIRS VNP43IA2 Quality (Collection 001)'
  });
});

Export.table.toDrive({
  collection: areaMeanQuality,
  description: 'HeardIsland_VIIRS_Quality_AreaMean',
  fileNamePrefix: 'heard_island_viirs_quality_area_mean',
  fileFormat: 'CSV',
  folder: driveFolder,
  selectors: ['date', 'quality_i1', 'quality_i2', 'quality_i3', 'valid_days_i1', 'valid_days_i2', 'valid_days_i3', 'snow_flag', 'land_water_type', 'product', 'source']
});

// Albedo area-mean export
var areaMeanAlbedo = viirsAlbedoProcessed.map(function(img) {
  var meanDict = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: glacierGeometry,
    scale: 375,
    maxPixels: 1e13,
    tileScale: 4
  });
  
  return ee.Feature(null, {
    date: ee.Date(img.get('system:time_start')).format('YYYY-MM-dd'),
    albedo_i1: meanDict.get('albedo_i1'),
    albedo_i2: meanDict.get('albedo_i2'),
    albedo_i3: meanDict.get('albedo_i3'),
    albedo_broadband: meanDict.get('albedo_broadband'),
    product: 'VNP43IA1',
    source: 'VIIRS VNP43IA1 BRDF (Collection 001)'
  });
});

Export.table.toDrive({
  collection: areaMeanAlbedo,
  description: 'HeardIsland_VIIRS_Albedo_AreaMean',
  fileNamePrefix: 'heard_island_viirs_albedo_area_mean',
  fileFormat: 'CSV',
  folder: driveFolder,
  selectors: ['date', 'albedo_i1', 'albedo_i2', 'albedo_i3', 'albedo_broadband', 'product', 'source']
});

// =============================================================================
// 6. YEARLY COMPOSITES
// =============================================================================
function createYearlyComposites(collection, productName, startYear, endYear) {
  var years = ee.List.sequence(startYear, endYear).getInfo();
  
  return ee.ImageCollection.fromImages(
    years.map(function(y) {
      var yearNum = ee.Number(y);
      var yearCol = collection.filter(ee.Filter.calendarRange(y, y, 'year'));
      var count = yearCol.size();
      
      return ee.Image(ee.Algorithms.If(
        count.gt(0),
        yearCol.mean()
          .set('year', yearNum)
          .set('product', productName)
          .set('system:time_start', ee.Date.fromYMD(yearNum, 1, 1).millis()),
        ee.Image.constant(0).rename('albedo').updateMask(0)
          .set('year', yearNum)
          .set('product', productName)
          .set('system:time_start', ee.Date.fromYMD(yearNum, 1, 1).millis())
      ));
    }).filter(function(img) {
      return img.select('albedo').bandNames().size().gt(0);
    })
  );
}

// Create yearly composites for VIIRS albedo
var yearlyAlbedo = createYearlyComposites(viirsAlbedoProcessed.select('albedo_broadband'), 'VNP43IA1', 2012, 2025);

var year2020Albedo = yearlyAlbedo.filter(ee.Filter.eq('year', 2020)).first();
Map.addLayer(year2020Albedo, albedoVis, 'VIIRS Albedo - 2020 Yearly Mean');

// =============================================================================
// 7. YEARLY RASTER EXPORTS
// =============================================================================
ee.List.sequence(2012, 2025).getInfo().forEach(function(y) {
  var img = ee.Image(yearlyAlbedo.filter(ee.Filter.eq('year', y)).first());
  Export.image.toDrive({
    image: img,
    description: 'HeardIsland_VIIRS_Albedo_Yearly_' + y,
    fileNamePrefix: 'heard_island_viirs_albedo_yearly_' + y,
    region: glacierGeometry,
    scale: 375,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    folder: driveFolder
  });
});

// =============================================================================
// 8. SUMMARY (console-safe)
// =============================================================================
safePrintIC('VNP43IA1 BRDF (raw)', viirsBRDF);
safePrintIC('VNP43IA2 Quality (raw)', viirsQuality);
// =============================================================================
// END OF SCRIPT
// =============================================================================
