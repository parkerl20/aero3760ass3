import ee
import geemap
import rasterio
import matplotlib.pyplot as plt
import imageio.v2 as imageio
from PIL import Image


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


def calculate_ndvi(image):
    # Calculates the NDVI index, which is a function of both the NIR and RED bands
    ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI')
    return image.addBands(ndvi)


def initialise_credentials():
    # Initialise credentials
    service_account = "spacey@spacey-401400.iam.gserviceaccount.com"
    key_path = "gee/key.json"    
    credentials = ee.ServiceAccountCredentials(service_account, key_path)
    ee.Initialize(credentials=credentials)


def dataset_to_tif(image, file_path, scale, band, rectangle_bounds):
    """
    Export an image to a TIF file with a specified scale and region of interest using the geemap.ee_export_image() function.

    Args:
    image (ee.Image): The Earth Engine image to be exported.
    file_path (str): The file path where the TIF file will be saved.
    scale (int): The scale at which to export the image.
    lat (float): Latitude coordinate of the region of interest.
    lon (float): Longitude coordinate of the region of interest.
    band (str): The band of the image to export.

    Returns:
    str: A message indicating the success or failure of the export process.
    """
    try:
        roi = ee.Geometry.Rectangle(rectangle_bounds) # Example region coordinates
        image = image.select(band)
        geemap.ee_export_image(image, filename=file_path, scale=scale, crs='EPSG:4326', region=roi, file_per_band=False)
        print(f"Image exported successfully to {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")
    


def tif_to_png(file_paths, i, gif_path):


    images = []

    # Create a figure and plot the image
    fig, ax = plt.subplots(figsize=(10, 10))

    for i, file_path in enumerate(file_paths):
        with rasterio.open(file_path) as src:
            img = src.read(2)  # Replace '1' with the band you want to visualize
            transform = src.transform
            extent = [transform[2], transform[2] + transform[0] * src.width,
                    transform[5] + transform[4] * src.height, transform[5]]

            ax.imshow(img, extent=extent, alpha=0.5)

            # Set axis labels
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')

            # Limits of NSW
            ax.set_xlim(141.0000, 156.6400)
            ax.set_ylim(-37.4000, -28.6000)

            # Show the plot
            print("Plot saved!")
            plt.savefig(f'../figures/dataset_GEE_export_{i+1}.png')

            images.append(imageio.imread(f'../figures/dataset_GEE_export_{i+1}.png'))

    # Save the GIF
    imageio.mimsave(gif_path, images, duration=100, loop=0)
    print(f"GIF saved at {gif_path}")


def generate_gif(lon_lat):

    initialise_credentials()

    dataset = (
        ee.ImageCollection('COPERNICUS/S2_SR')
        .filterDate("2019-12-01", "2020-01-31")
        # Pre-filter to get less cloudy granules.
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
        .map(mask_s2_clouds)
        .map(calculate_ndvi) # Calculates the NDVI index
    )

    file_paths = []
    gif_path = '../figures/dataset_spatiotemporal_swath.gif'

    for i in range(0, len(lon_lat)):

        rectangle_bounds = [lon_lat[i][0] - 3.1, lon_lat[i][1] + 3.1, lon_lat[i][0] + 3.1, lon_lat[i][1] - 3.1] # Scale 800

        coverage = ee.Geometry.Rectangle(rectangle_bounds)

        infrared_clipped = dataset.mean().clip(coverage)

        # Specify the file path, scale, and region for export
        file_path = f'tifs/image-{i+1}.tif'
        scale = 800
        band = ['B8', 'B4', 'B3'] # Bands for infrared

        # tif --> png --> saved in files
        dataset_to_tif(infrared_clipped, file_path, scale, band, rectangle_bounds)
        file_paths.append(file_path)

    tif_to_png(file_paths, i, gif_path)

    return 0

