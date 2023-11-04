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
    key_path = "key.json"    
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

            # Add a colorbar
            # cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            # cbar.set_label('Value')

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


def st_main():

    lon_lat = [[148.6174199280952, -28.83543367853245], [149.92298543700326, -29.27253474647698], 
               [151.23867040151083, -29.694669800958295], [152.56431370584312, -30.10151175342973], 
               [153.89972151100193, -30.49274096785361], [155.24466667775175, -30.868046182909563], 
               [156.5988883533659, -31.227125459097678], [142.00831136929904, -33.21417299045468], 
               [143.4220963939283, -33.43790770103715], [144.8414445346663, -33.643314257897444], 
               [146.26578656207727, -33.83023255351011], [147.69452894296646, -33.9985216486313], 
               [149.1270561076186, -34.14806043632942], [150.562732924162, -34.27874822682028], 
               [152.0009073605424, -34.39050524578126], [153.44091331120535, -34.483273039666834], 
               [154.88207356253966, -34.557014782531574], [156.32370286851835, -34.61171547997106], 
               [142.06694871160516, -34.482792569382866], [143.4888471513621, -34.39113317793818], 
               [144.90533668993, -34.281222600031846], [146.31581748390514, -34.153248331004754], 
               [147.71971323569647, -34.007415349490294], [149.1164732503732, -33.843945211655964], 
               [150.50557425767715, -33.66307509206819], [151.88652199062489, -33.465056780641845], 
               [153.25885251589696, -33.25015564524047], [154.62213331486942, -33.01864956944403], 
               [155.9759641175872, -32.770827874806066], [142.16001498019, -30.317887934379836], 
               [143.41869990068852, -29.940425851962235], [144.66574542546567, -29.549864661282708], 
               [145.90109457113422, -29.146551736084454], [147.12471770123528, -28.73083458991063]]

    # lon_lat = [[142.80918267990074, -28.900726724884922], [144.1162472170203, -29.335644741093237], [145.43340951232395, -29.75554814475814], [146.76050367340613, -30.160110887976725], [148.09733099005928, -30.549014509388723], [149.44365938028915, -30.921949060387036], [150.79922300261782, -31.27861405541453], [152.16372205140075, -31.618719440453113], [153.53682275083318, -31.941986573130496], [154.91815756187881, -32.24814920723015], [156.30732561451447, -32.53695447380844], [141.90540579012327, -34.02172220762589], [143.33843707707317, -34.16849087667674], [144.7745229944953, -34.29639547718888], [146.21300925680197, -34.40535935133647], [147.65322793353445, -34.49532721355476], [149.09450043288436, -34.56626544246081], [150.53614059101346, -34.61816226934489], [151.97745783622793, -34.65102786022855], [153.41776039560753, -34.66489428976839], [154.85635851085505, -34.65981540661138], [156.29256762994774, -34.63586659114664], [141.92600655880432, -33.98444505825686], [143.32167090439737, -33.818397396303084], [144.70960160257962, -33.63498593737796], [146.08930877567596, -33.43446453668128], [147.46033305450325, -33.21710047227287], [148.82224662279623, -32.98317337997729], [150.17465401519112, -32.73297417412627], [151.51719267469443, -32.46680396307726], [152.84953327837428, -32.1849729679899], [154.17137984246347, -31.887799452791562], [155.48246962012402, -31.5756086726344], [141.3040634016958, -28.66859050326423]]


    print("LON_LAT:", lon_lat)

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
        # rectangle_bounds = [141.115256 + 2*i, -28.356159 - 1*i, 147.315256 + 2*i, -34.956159 - 1*i] # Scale 800

        # rectangle_bounds = [151.115256, -33.756159, 151.315256, -33.956159] # Centered at the Opera House, scale 60
        # rectangle_bounds = [151.115256, -33.756159, 151.125256, -33.766159] # Very very pixellated
        # rectangle_bounds = [150.115256, -32.756159, 152.315256, -34.956159] # Centered at the Opera House, scale 200
        # rectangle_bounds = [141.0000, -29.1770, 153.6372, -37.5050] # Full NSW, supposedly. 

        # scale = (rectangle_bounds[2] - rectangle_bounds[0]) * 200

        # squares = ee.Geometry.Point([rectangle_bounds[0], rectangle_bounds[1]]).buffer(345088).bounds()
        # coverage = ee.Geometry.MultiPolygon([squares.coordinates()])

        coverage = ee.Geometry.Rectangle(rectangle_bounds)

        infrared_clipped = dataset.mean().clip(coverage)

        # Specify the file path, scale, and region for export
        file_path = f'tifs/image-{i+1}.tif'
        scale = 800
        '''
        I have found that scale = 60 is the smallest it goes. 
        Larger the scale, the more zoomed in, the smaller bits are downloaded.
        Smaller the scale, the more zoomed out, the more bits are downloaded.

        The scale represents the GSD. If scale = 20, the image has GSD of 20m
        '''
        band = ['B8', 'B4', 'B3'] # Bands for infrared

        # tif --> png --> saved in files
        dataset_to_tif(infrared_clipped, file_path, scale, band, rectangle_bounds)
        file_paths.append(file_path)

    tif_to_png(file_paths, i, gif_path)

    return 0


if __name__ == "__main__":
    st_main()
