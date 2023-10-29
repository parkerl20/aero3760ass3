import ee
import geemap
import rasterio
import matplotlib.pyplot as plt
import requests
from ee import batch
import os


def initialise_credentials():
    # Initialise credentials
    service_account = "spacey@spacey-401400.iam.gserviceaccount.com"
    key_path = "key.json"    
    credentials = ee.ServiceAccountCredentials(service_account, key_path)
    ee.Initialize(credentials=credentials)


def export_image_to_tif_with_region(image, file_path, scale, lat, lon, band):
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
        roi = ee.Geometry.Point(lon, lat).buffer(10000)  # Adjust buffer size as needed
        image = image.select(band)
        geemap.ee_export_image(image, filename=file_path, scale=scale, region=roi)
        return f"Image exported successfully to {file_path}"
    except Exception as e:
        return f"An error occurred: {e}"
    


def save_to_png(file_path):
    with rasterio.open(file_path) as src:
        # Read the image as a numpy array
        img = src.read(1)  # Replace '1' with the band you want to visualize

        # Get the spatial transform (geotransform)
        transform = src.transform

        # Get the spatial extent (bounding box)
        extent = [transform[2], transform[2] + transform[0] * src.width,
                transform[5] + transform[4] * src.height, transform[5]]

        # Create a figure and plot the image
        fig, ax = plt.subplots(figsize=(10, 10))
        im = ax.imshow(img, extent=extent)

        # Add a colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Value')

        # Set axis labels
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')

        # Show the plot
        print("Plot saved!")
        plt.show()
        plt.savefig('figures/sydney.png')


def main_dataset2image():

    # Define the region of interest (AOI) as a geometry
    aoi = ee.Geometry.Polygon(
        [[[-151.5, 33.5],
        [-149.4, 33.5],
        [-149.4, 35.6],
        [-151.5, 35.6]]])

    # Define the date range
    start_date = '2021-01-01'
    end_date = '2021-12-31'

    # Load the Landsat 8 dataset for the specified time period and location
    dataset = ee.ImageCollection('LANDSAT/LC08/C01/T1') \
        .filterBounds(aoi) \
        .filterDate(ee.Date(start_date), ee.Date(end_date)) \
        .sort('CLOUD_COVER') \
        .first()  # Get the image with the least cloud cover

    # Apply the clipToBoundsAndScale function with the input parameter
    clipped_image = dataset.geometry().clipToBoundsAndScale(aoi)

    # Define export parameters
    export_params = {
        'image': clipped_image,
        'description': 'Landsat8_Image',
        'scale': 30,  # Spatial resolution in meters per pixel
        'region': aoi.getInfo()['coordinates'],
        'fileFormat': 'GeoTIFF'  # Export format (e.g., GeoTIFF)
    }

    # # Load the Landsat 8 dataset for the specified time period and location
    # dataset = ee.ImageCollection('LANDSAT/LC08/C01/T1') \
    #     .filterBounds(aoi) \
    #     .filterDate(ee.Date(start_date), ee.Date(end_date)) \
    #     .sort('CLOUD_COVER') \
    #     .first()  # Get the image with the least cloud cover

    # # Select the bands you want to export (e.g., Red, Green, Blue)
    # bands_to_export = ['B4', 'B3', 'B2']  # RGB bands

    # # Define export parameters
    # export_params = {
    #     'image': dataset.select(bands_to_export),
    #     'description': 'Landsat8_Image',
    #     'scale': 30,  # Spatial resolution in meters per pixel
    #     'region': aoi.getInfo()['coordinates'],
    #     'fileFormat': 'GeoTIFF'  # Export format (e.g., GeoTIFF)
    # }

    # Initiate the export task
    task = batch.Export.image.toDrive(**export_params)

    # Monitor the task status
    print('Exporting...', task)

    # Optionally, wait for the task to complete
    while task.active():
        pass

    # Check the task status
    print('Task Status:', task.status())

    # Get the download URL for the exported image
    download_url = dataset.getDownloadURL({
        'name': 'Landsat8_Image',
        'scale': 30,
        'region': aoi
    })

    # Define the local output path for the downloaded image
    output_path = 'tifs/Landsat8_Image.tif'

    # Download the exported image from the URL
    response = requests.get(download_url)
    with open(output_path, 'wb') as f:
        f.write(response.content)

    # Check if the image was successfully downloaded
    if os.path.exists(output_path):
        print('Image downloaded successfully.')
    else:
        print('Image download failed.')

    return 0



def save_to_tiff(dataset, aoi, output_path):
    # Clip the dataset to the region of interest
    clipped_image = dataset.map(lambda image: image.clip(aoi))

    # Get the download URL for the clipped image
    download_url = clipped_image.select('NDVI').getDownloadURL({
        'name': 'NDVI_Image',
        'scale': 10,  # Spatial resolution in meters per pixel
        'region': aoi
    })

    # Download the image from the URL
    response = requests.get(download_url)

    # Save the image to the specified output path
    with open(output_path, 'wb') as f:
        f.write(response.content)

    # Check if the image was successfully downloaded
    if os.path.exists(output_path):
        print(f'Image downloaded successfully at {output_path}.')
    else:
        print('Image download failed.')


    return 0





def main():

    initialise_credentials()

    # Usage
    start_date = '2021-01-01'
    end_date = '2021-12-31'
    aoi = ee.Geometry.Polygon(
            [[[-122.45, 37.74],
            [-122.4, 37.74],
            [-122.4, 37.8],
            [-122.45, 37.8]]])  # Example region of interest

    # Define the dataset
    dataset = (
        ee.ImageCollection('COPERNICUS/S2_SR')
        .filterDate(start_date, end_date)
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
        # Add your custom functions here, for example:
        # .map(mask_s2_clouds)
        # .map(calculate_ndvi)
    )

    # Specify the output path for the TIFF file
    output_path = 'tifs/output.tif'

    # Call the function to save the dataset to a TIFF file
    save_to_tiff(dataset, aoi, output_path)

    # main_dataset2image()

    # # Example usage
    # # Load an image (you should replace this with your own Earth Engine image)
    # image = ee.Image('USGS/SRTMGL1_003')

    # # Specify the file path, scale, and region for export
    # file_path = 'tifs/image.tif'
    # scale = 30
    # lat = -33.7749  # Example latitude
    # lon = 151.4194  # Example longitude
    # band = 'elevation'  # Example band name

    # # Export the image
    # export_result = export_image_to_tif_with_region(image, file_path, scale, lat, lon, band)
    # print(export_result)

    # save_to_png(file_path)

    return 0


if __name__ == "__main__":
    main()
