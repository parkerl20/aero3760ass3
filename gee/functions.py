import ee
import os
import webbrowser
import geemap
import csv

def initialise_credentials():
    # Initialise credentials
    service_account = "spacey@spacey-401400.iam.gserviceaccount.com"
    key_path = "gee/key.json"    
    credentials = ee.ServiceAccountCredentials(service_account, key_path)
    ee.Initialize(credentials=credentials)


def create_map(map_object, map_name):
    # Directory
    directory = os.path.join("gee", "maps")
    file_extension = ".html"
    
    # Build the full file path
    file_path = os.path.join(directory, map_name + file_extension)

    # Save the map as a html
    map_object.save(file_path)

def show_map(map_name):
    # Directory
    directory = os.path.join("gee", "maps")
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


def S2A_NDVI(start_date: str, end_date: str):
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
        .map(calculate_ndvi) # Calculates the NDVI index
    )

    # Selecting NDVI
    mean_ndvi = dataset.mean().select('NDVI')

    # Visualisation parameter
    ndvi_vis = {
        'min': -1,
        'max': 1,
        'palette': ['blue', 'white', 'green']
    }

    # Map initialisation
    Map = geemap.Map() 
    Map.set_center(146.9211, -31.2532, 6) # Center of nsw
    Map.add_ee_layer(mean_ndvi, ndvi_vis, 'NDVI')

    # Locations to take NDVI data
    locations = ee.FeatureCollection([
        ee.Feature(ee.Geometry.Point([146.9211, -31.2532]), {'name': 'NSW Centre'}),
        ee.Feature(ee.Geometry.Point([151.2093, -33.8688]), {'name': 'Sydney'}),
        ee.Feature(ee.Geometry.Point([152.1895, -30.8487]), {'name': 'Carrai Creek'}),
        ee.Feature(ee.Geometry.Point([151.2642, -33.6184]), {'name': 'Ku ring gai chase park'})
    ])

    # Get NDVI values 
    ndvi_values = extract_ndvi_values(mean_ndvi, locations)
    print(ndvi_values)

    # Write the NDVI values to a csv
    write_to_csv(ndvi_values, 'ndvi_values.csv', 'NDVI Index')

    return Map


def S2A_coverage(start_date: str, end_date: str):
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
        .map(calculate_ndvi) # Calculates the NDVI index
    )

    # Selecting NDVI
    mean_ndvi = dataset.mean().select('NDVI')

    # Visualisation parameter
    ndvi_vis = {
        'min': -1,
        'max': 1,
        'palette': ['blue', 'white', 'green']
    }

    # Map initialisation
    Map = geemap.Map() 
    Map.set_center(146.9211, -31.2532, 6) # Center of nsw

    # Coverage region
    points_of_interest = [
        [151.2093, -33.8688],  # Sydney
        [150.8931, -34.4278],  # Wollongong
        [151.7817, -32.9275]   # Newcastle
    ]

    # Coverage points
    multipoint = ee.Geometry.MultiPoint(points_of_interest)
    coverage = multipoint.buffer(2000)
    mean_ndvi_clipped = mean_ndvi.clip(coverage)
    Map.add_ee_layer(mean_ndvi_clipped, ndvi_vis, "NDVI")

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


def calculate_ndvi(image):
    # Calculates the NDVI index, which is a function of both the NIR and RED bands
    ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI')
    return image.addBands(ndvi)


def dynamic_world():
    # Construct a collection of corresponding Dynamic World and Sentinel-2 for
    # inspection. Filter the DW and S2 collections by region and date.
    START = ee.Date('2021-04-02')
    END = START.advance(1, 'day')

    # Get the geometry of New South Wales
    nsw = ee.FeatureCollection("FAO/GAUL/2015/level1").filter(
            ee.Filter.eq('ADM1_NAME', 'New South Wales')
    )

    col_filter = ee.Filter.And(
        ee.Filter.bounds(nsw),
        ee.Filter.date(START, END),
    )

    dw_col = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1').filter(col_filter)
    s2_col = ee.ImageCollection('COPERNICUS/S2').filter(col_filter)

    # Join corresponding DW and S2 images (by system:index).
    dw_s2_col = ee.Join.saveFirst('s2_img').apply(
        dw_col,
        s2_col,
        ee.Filter.equals(leftField='system:index', rightField='system:index'),
    )

    # Extract an example DW image and its source S2 image.
    dw_image = ee.Image(dw_s2_col.first())
    s2_image = ee.Image(dw_image.get('s2_img'))

    # Create a visualization that blends DW class label with probability.
    # Define list pairs of DW LULC label and color.
    CLASS_NAMES = [
        'water',
        'trees',
        'grass',
        'flooded_vegetation',
        'crops',
        'shrub_and_scrub',
        'built',
        'bare',
        'snow_and_ice',
    ]

    VIS_PALETTE = [
        '419bdf',
        '397d49',
        '88b053',
        '7a87c6',
        'e49635',
        'dfc35a',
        'c4281b',
        'a59b8f',
        'b39fe1',
    ]

    # Create an RGB image of the label (most likely class) on [0, 1].
    dw_rgb = (
        dw_image.select('label')
        .visualize(min=0, max=8, palette=VIS_PALETTE)
        .divide(255)
    )

    # Get the most likely class probability.
    top1_prob = dw_image.select(CLASS_NAMES).reduce(ee.Reducer.max())

    # Create a hillshade of the most likely class probability on [0, 1]
    top1_prob_hillshade = ee.Terrain.hillshade(top1_prob.multiply(100)).divide(255)

    # Combine the RGB image with the hillshade.
    dw_rgb_hillshade = dw_rgb.multiply(top1_prob_hillshade)

    # Display the Dynamic World visualization with the source Sentinel-2 image.
    Map = geemap.Map()
    Map.set_center(20.6729, 52.4305, 12)
    Map.add_layer(
        s2_image,
        {'min': 0, 'max': 3000, 'bands': ['B4', 'B3', 'B2']},
        'Sentinel-2 L1C',
    )
    Map.add_layer(
        dw_rgb_hillshade,
        {'min': 0, 'max': 0.65},
        'Dynamic World V1 - label hillshade',
    )

    return Map


def fires():
    # Get the geometry of New South Wales
    nsw = ee.FeatureCollection("FAO/GAUL/2015/level1").filter(
            ee.Filter.eq('ADM1_NAME', 'New South Wales')
    )

    dataset = (
        ee.ImageCollection('FIRMS')
        .filterBounds(nsw)
        .filterDate('2019-12-01', '2020-02-01')
    )
    fires = dataset.select('T21')
    firesVis = {
        "min": 325.0,
        "max": 400.0,
        "palette": ['red', 'orange', 'yellow'],
    }

    # Map initialisation
    Map = geemap.Map() 
    Map.set_center(146.9211, -31.2532, 6) # Center of nsw
    Map.addLayer(fires, firesVis, 'Fires')

    return Map


def extract_ndvi_values(image, feature_collection):
    # Get a feature collection from the imnage
    sampled_points = image.select(['NDVI']).sampleRegions(
        collection=feature_collection, 
        properties=['name'], 
        scale=10  
    )
    
    # Convert the feature collection to a list and evaluate it.
    info = sampled_points.getInfo()['features']
    
    # Extract NDVI values from the returned info and store them in a dictionary.
    ndvi_values = {}
    for feat in info:
        ndvi_values[feat['properties']['name']] = feat['properties']['NDVI']
        
    return ndvi_values


def write_to_csv(data, filename, band):
    # Writes band data to a csv
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Location", band])  # Writing headers
        for location, value in data.items():
            writer.writerow([location, value])