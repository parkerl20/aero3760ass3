import ee
import os
import webbrowser
import geemap
import csv
from spacesim import orbital_transforms as ot

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


def S2A_coverage(start_date: str, end_date: str, lon_lat, circle_radius):
   # Sentinel-2A satellite
    dataset = (
        ee.ImageCollection('COPERNICUS/S2_SR')
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

    # Red for alternate plotting
    infrared_vis = {
        'min': 0.0,
        'max': 0.3,
        'bands': ['B8', 'B4', 'B3']
    }

    # Map initialisation
    Map = geemap.Map() 
    Map.set_center(146.9211, -31.2532, 6) # Center of nsw

    # Coverage points
    # multipoint = ee.Geometry.MultiPoint(lon_lat)

    point = ee.Geometry.Point(lon_lat[6]) # Point near sydney
    # coverage = point.buffer(circle_radius)

    # squares = [ee.Geometry.Point(coords).buffer(circle_radius).bounds() for coords in lon_lat]
    # coverage = ee.Geometry.MultiPolygon([square.coordinates() for square in squares])

    squares = ee.Geometry.Point(lon_lat[6]).buffer(circle_radius).bounds()
    coverage = ee.Geometry.MultiPolygon([squares.coordinates()])

    # coverage = multipoint.buffer(circle_radius)
    # mean_ndvi_clipped = mean_ndvi.clip(coverage)
    infrared_clipped = dataset.mean().clip(coverage)
    # Map.add_ee_layer(mean_ndvi_clipped, ndvi_vis, "NDVI")
    Map.add_ee_layer(infrared_clipped, infrared_vis, "Coverage")

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



def plot_red_points(lon_lat, circle_radius):
    # Red color visualization
    red_vis = {
        'palette': ['red']
    }
    # Adding the opera house as a point
    lon_lat.append([151.2153, -33.8568])
    
    # Create a MultiPoint geometry from the coordinates
    multipoint = ee.Geometry.MultiPoint(lon_lat)
    
    # Paint the multipoint geometry on an image
    buffered_multipoint = multipoint.buffer(circle_radius)  # Radius of circle
    red_points_image = ee.Image().paint(buffered_multipoint, 0)

    # Map initialization (empty map)
    Map = geemap.Map()
    Map.set_center(146.9211, -31.2532, 6)  # Center of NSW
    Map.add_ee_layer(red_points_image, red_vis, 'Red Points')

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

def eci_to_llh(r_eci, t_eci, epoch):
    # 2D array for latitude longitude storage
    lon_lat = [[0.0, 0.0] for i in range(len(t_eci))]

    # Convert eci to lat long
    for i in range(len(t_eci)):
        r_ecef = ot.ECI_to_ECEF(r_eci[:, i], t_eci[i], epoch)
        lat_i, lng_i, _ = ot.ECEF_to_LLH(r_ecef).flatten()
        lon_lat[i] = [lng_i, lat_i]

    return lon_lat


def eci_to_llh_nsw(r_eci, t_eci, epoch, num_points):
    # 2D array for latitude longitude storage
    lon_lat = [[0.0, 0.0] for i in range(len(t_eci))]

    # Convert eci to lat long
    for i in range(len(t_eci)):
        r_ecef = ot.ECI_to_ECEF(r_eci[:, i], t_eci[i], epoch)
        lat_i, lng_i, _ = ot.ECEF_to_LLH(r_ecef).flatten()
        lon_lat[i] = [lng_i, lat_i]

    # Bounds of NSW
    nsw_bounds = {
        "lat_min": -37.4,
        "lat_max": -28.6,
        "lon_min": 141,
        "lon_max": 156.64
    }

    # Only consider points in NSW
    filtered_lon_lat = [point for point in lon_lat 
                        if nsw_bounds["lat_min"] <= point[1] <= nsw_bounds["lat_max"] 
                        and nsw_bounds["lon_min"] <= point[0] <= nsw_bounds["lon_max"]]
    
    # When extra points are not needed for visualisation
    if(num_points == 0):
        return filtered_lon_lat


    # Adding interpolated points between consecutive points
    smooth_lon_lat = []
    for i in range(len(filtered_lon_lat) - 1):
        smooth_lon_lat.append(filtered_lon_lat[i])
        smooth_lon_lat.extend(interpolate_points(filtered_lon_lat[i], filtered_lon_lat[i+1], num_points)) 
    if filtered_lon_lat:
        smooth_lon_lat.append(filtered_lon_lat[-1])

    return smooth_lon_lat



def interpolate_points(A, B, num_points):
    # Add in points between lon_lat for a smoother curve    
    x1, y1 = A
    x2, y2 = B
    
    # Delta
    x_delta = (x2 - x1) / (num_points + 1)
    y_delta = (y2 - y1) / (num_points + 1)
    
    # Generate interpolated points
    interpolated_points = []
    for i in range(1, num_points + 1):
        new_x = x1 + i * x_delta
        new_y = y1 + i * y_delta
        interpolated_points.append((new_x, new_y))
        
    return interpolated_points
