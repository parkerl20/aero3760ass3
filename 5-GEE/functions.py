import ee
import os
import webbrowser
import geemap

def initialise_credentials():
    # Initialise credentials
    service_account = "spacey@spacey-401400.iam.gserviceaccount.com"
    key_path = "key.json"    
    credentials = ee.ServiceAccountCredentials(service_account, key_path)
    ee.Initialize(credentials=credentials)


def create_map(map_object, map_name):
    # Directory
    directory = os.path.join("maps")
    file_extension = ".html"
    
    # Build the full file path
    file_path = os.path.join(directory, map_name + file_extension)

    # Save the map as a html
    map_object.save(file_path)

def show_map(map_name):
    # Directory
    directory = os.path.join("maps")
    file_extension = ".html"
    
    # Build the full file path
    file_path = os.path.join(directory, map_name + file_extension)

    webbrowser.open('file:///' + os.path.realpath(file_path))

def mask_s2_clouds(image):
  """Masks clouds in a Sentinel-2 image using the QA band.

  Args:
      image (ee.Image): A Sentinel-2 image.

  Returns:
      ee.Image: A cloud-masked Sentinel-2 image.
  """
  qa = image.select('QA60')

  # Bits 10 and 11 are clouds and cirrus, respectively.
  cloud_bit_mask = 1 << 10
  cirrus_bit_mask = 1 << 11

  # Both flags should be set to zero, indicating clear conditions.
  mask = (
      qa.bitwiseAnd(cloud_bit_mask)
      .eq(0)
      .And(qa.bitwiseAnd(cirrus_bit_mask).eq(0))
  )

  return image.updateMask(mask).divide(10000)


def initial_map():
    # Set some points of interest
    # Mackay, Qld
    lat = -21.1434
    lon = 149.1868

    # Get URL for our image
    Data_Set_Clearing = ee.ImageCollection("WRI/GFW/FORMA/raw_output_firms")

    # Create an interactive map
    # Map = folium.Map(center=[lat, lon], zoom_start=7)
    Map = geemap.Map(center=[lat, lon], zoom=7)

    # Define the dataset
    dataset = ee.ImageCollection('WRI/GFW/FORMA/raw_output_firms')
    dataset.filter(ee.Filter.date('2018-08-01', '2018-08-15'))

    # Select the 'nday' band
    percentageOfClearing = dataset.select('nday')

    # Visualization parameters
    visParams = {
    'min': 0.0,
    'max': 0.01
    }

    # Add layer to map
    Map.addLayer(percentageOfClearing, visParams, 'Percentage of clearing')

    return Map

def S2A(start_date: str, end_date: str):
    # Get the geometry of New South Wales
    nsw = ee.FeatureCollection("FAO/GAUL/2015/level1").filter(
            ee.Filter.eq('ADM1_NAME', 'New South Wales')
    )
 
   # Sentinel-2A satellite
    dataset = (
        ee.ImageCollection('COPERNICUS/S2_SR')
        .filterBounds(nsw)
        .filterDate(start_date, end_date)
        # Pre-filter to get less cloudy granules.
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
        .map(mask_s2_clouds)
    )

    # Visualisation parameter
    visualization = {
        'min': 0.0,
        'max': 0.3,
        'bands': ['B8', 'B4', 'B3'],
    }

    # Map initialisation
    Map = geemap.Map() 
    Map.set_center(146.9211, -31.2532, 6) # Center of nsw
    Map.add_ee_layer(dataset.mean(), visualization, 'Infrared')

    return Map


def elevation_5m():
    # Get the geometry of New South Wales
    nsw = ee.FeatureCollection("FAO/GAUL/2015/level1").filter(
            ee.Filter.eq('ADM1_NAME', 'New South Wales')
    )

    # 5m elevation dataset
    dataset = (
        ee.ImageCollection('AU/GA/AUSTRALIA_5M_DEM')
        .filterBounds(nsw)
    )

    # Visualisation parameter
    elevation = dataset.select('elevation')
    elevation_visualization = {
        'min': 0.0,
        'max': 150.0,
        'palette': ['0000ff', '00ffff', 'ffff00', 'ff0000', 'ffffff'],
    }

    # Map initialisation
    Map = geemap.Map() 
    Map.addLayer(elevation, elevation_visualization, 'Elevation')

    return Map


def surface_temperature():
    # Get the geometry of New South Wales
    nsw = ee.FeatureCollection("FAO/GAUL/2015/level1").filter(
            ee.Filter.eq('ADM1_NAME', 'New South Wales')
    )

    # 5m elevation dataset
    dataset = (
        ee.ImageCollection('JAXA/GCOM-C/L3/LAND/LST/V3')
        .filterBounds(nsw)
        .filterDate('2019-12-01', '2020-02-01')
        .filter(ee.Filter.eq('SATELLITE_DIRECTION', 'D'))
    )

    # Multiply with slope coefficient
    dataset = dataset.mean().multiply(0.02)

    # Visualisation parameter
    visualization = {
        'bands': ['LST_AVE'],
        'min': 250,
        'max': 316,
        'palette': ['040274','040281','0502a3','0502b8','0502ce','0502e6',
                    '0602ff','235cb1','307ef3','269db1','30c8e2','32d3ef',
                    '3be285','3ff38f','86e26f','3ae237','b5e22e','d6e21f',
                    'fff705','ffd611','ffb613','ff8b13','ff6e08','ff500d',
                    'ff0000','de0101','c21301','a71001','911003'],
    }

    # Map initialisation
    Map = geemap.Map() 
    Map.addLayer(dataset, visualization, 'Elevation')

    return Map
