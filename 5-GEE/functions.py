import ee
import os
import webbrowser
import geemap

def initialise_credentials():
    # Initialise credentials
    service_account = "spacey@spacey-401400.iam.gserviceaccount.com"
    key_path = "5-GEE/key.json"    
    credentials = ee.ServiceAccountCredentials(service_account, key_path)
    ee.Initialize(credentials=credentials)


def create_map(map_object, map_name):
    # Directory
    directory = os.path.join("5-GEE", "maps")
    file_extension = ".html"
    
    # Build the full file path
    file_path = os.path.join(directory, map_name + file_extension)

    # Save the map as a html
    map_object.save(file_path)

def show_map(map_name):
    # Directory
    directory = os.path.join("5-GEE", "maps")
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

def S2A():
   # Sentinel-2A satellite
    dataset = (
        ee.ImageCollection('COPERNICUS/S2_SR')
        .filterDate('2020-01-01', '2020-01-30')
        # Pre-filter to get less cloudy granules.
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
        .map(mask_s2_clouds)
    )

    # Visualisation parameter
    visualization = {
        'min': 0.0,
        'max': 0.3,
        # 'bands': ['B8', 'B4', 'B3'],
        'bands': ['B4', 'B3', 'B2'],
    }

    # Map initialisation
    Map = geemap.Map() 
    Map.set_center(146.9211, -31.2532, 6) # Center of nsw
    Map.add_ee_layer(dataset.mean(), visualization, 'Infrared')

    return Map