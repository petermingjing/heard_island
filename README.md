# Data Directory

This directory contains essential data files for the Heard Island Glacier Albedo analysis.

## Directory Structure

```
data/
├── gee_export/
│   ├── HIG_VIIRS_Albedo/
│   │   └── heard_island_viirs_albedo_area_mean.csv  # VIIRS time series
│   └── heard_island_era5_land_climate_daily.csv      # Climate data
└── heard_island_srtm_dem_resampling_viirs.tif        # DEM
```

## Essential Data Files

### 1. VIIRS Albedo Time Series

**File:** `gee_export/HIG_VIIRS_Albedo/heard_island_viirs_albedo_area_mean.csv`

- **Description:** Glacier-wide area-mean albedo time series
- **Format:** CSV
- **Columns:**
  - `date`: Observation date (YYYY-MM-DD)
  - `albedo`: White-sky albedo (unitless, 0-1)
- **Temporal Coverage:** January 18, 2012 to May 31, 2024
- **Number of Observations:** 4,466
- **Source:** NASA VIIRS VNP43 product
- **Size:** ~590 KB

### 2. Climate Data

**File:** `gee_export/heard_island_era5_land_climate_daily.csv`

- **Description:** Daily climate reanalysis data
- **Format:** CSV
- **Columns:**
  - `date`: Date (YYYY-MM-DD)
  - `temperature_2m`: 2-meter air temperature (°C)
  - `precipitation`: Total precipitation (mm)
  - `snowfall`: Snowfall water equivalent (mm)
  - `surface_solar_radiation_downwards`: Surface solar radiation (W/m²)
- **Temporal Coverage:** January 1, 2012 to May 31, 2024
- **Source:** ECMWF ERA5-Land reanalysis
- **Size:** ~500 KB

### 3. Digital Elevation Model

**File:** `heard_island_srtm_dem_resampling_viirs.tif`

- **Description:** SRTM DEM resampled to match VIIRS resolution
- **Format:** GeoTIFF
- **Spatial Resolution:** 375-500 m (matches VIIRS)
- **Projection:** EPSG:4326 (WGS84)
- **Units:** Meters above sea level
- **Source:** SRTM (Shuttle Radar Topography Mission)
- **Size:** ~26 KB

## Regenerating Data Files

### Annual TIF Rasters

Annual albedo rasters (`heard_island_viirs_albedo_yearly_YYYY.tif`) are not included in this repository but can be regenerated:

1. Open Google Earth Engine Code Editor
2. Copy script from `scripts/viirs.js`
3. Run script to export annual rasters
4. Download from Google Drive
5. Place in `data/gee_export/HIG_VIIRS_Albedo/`

### Climate Data

To regenerate climate data:

1. Open Google Earth Engine Code Editor
2. Copy script from `scripts/climate.js`
3. Run script to export climate data
4. Download from Google Drive
5. Place in `data/gee_export/`

### DEM Data

To regenerate DEM:

1. Open Google Earth Engine Code Editor
2. Copy script from `scripts/srtm_export.js`
3. Run script to export DEM
4. Download from Google Drive
5. Place in `data/`

## Data Sources

### Original Data Access

- **VIIRS Data:** Google Earth Engine (VNP43 product)
  - Collection: `NOAA/VIIRS/001/VNP43IA1` (BRDF parameters)
  - Collection: `NOAA/VIIRS/001/VNP43IA2` (Quality flags)

- **ERA5-Land Data:** Copernicus Climate Data Store
  - Product: ERA5-Land
  - Access: https://cds.climate.copernicus.eu/

- **SRTM DEM:** USGS EarthExplorer
  - Product: SRTM 1 Arc-Second
  - Access: https://earthexplorer.usgs.gov/

## Data Quality

### VIIRS Data Quality
- **Temporal Availability:** 98.9% (excellent)
- **Spatial Coverage:** 44.04% (moderate, typical for polar regions)
- **Quality Flags:** BRDF quality flags used for filtering
- **Validation:** Uses validated BRDF/albedo parameters

### Climate Data Quality
- **Source:** ERA5-Land reanalysis (state-of-the-art)
- **Spatial Resolution:** ~9 km (interpolated to point location)
- **Temporal Completeness:** 100% (no missing days)
- **Validation:** Validated against observations globally

## Excluded Data

The following data are excluded from this repository:

- **Archive directory** (`gee_export/archive/`): Contains excluded sensor data (MODIS, Landsat, Sentinel-2)
- **Annual TIF rasters**: Can be regenerated from GEE scripts
- **Original DEM**: Only resampled version included

See `.gitignore` for complete list of excluded files.

## Data Citation

When using this data, please cite:

1. **Original Data Sources:**
   - VIIRS: NASA VIIRS VNP43 product
   - ERA5-Land: Muñoz-Sabater et al. (2021)
   - SRTM: USGS SRTM

2. **This Dataset:**
   - Ming, J. (2024). Heard Island Glacier Albedo Dataset (2012-2024). GitHub. https://github.com/petermingjing/HIG

## Additional Information

For detailed data descriptions, see `docs/DATA_DESCRIPTION.md`.

